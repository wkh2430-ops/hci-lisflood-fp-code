"""Deterministic passive advection diagnostic, not calibrated tracer transport.
1000 equal-weight particles per cohort, positive-flux-weighted release locations.
Midpoint RK2, bilinear spatial and linear temporal interpolation, <=10 s steps.
All fractions use ALL released particles. Numerical exits are reported separately.
"""
from pathlib import Path
import json, argparse, time, os, sys
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent; O=P/'analysis'
parser=argparse.ArgumentParser(); parser.add_argument('--step',type=float,default=10); parser.add_argument('--cohort-limit',type=int); args=parser.parse_args()
default_run=args.step==10 and args.cohort_limit is None
if default_run:
    if (O/'particle_complete.txt').exists():
        print('Completed default particle result already exists; reusing it.',flush=True); sys.exit(0)
    try:
        with (O/'particle_running.lock').open('x') as f: f.write(str(os.getpid()))
    except FileExistsError:
        deadline=time.monotonic()+1800
        while not (O/'particle_complete.txt').exists():
            if time.monotonic()>deadline: raise RuntimeError('Existing particle run has not completed; inspect lock and run log')
            time.sleep(2)
        print('Concurrent default particle run completed.',flush=True); sys.exit(0)
g=json.loads((P/'experiment_manifest.json').read_text(encoding='utf-8'))['grid']
dx=g['dx']; left=g['left']; top=g['bottom']+g['ny']*dx
H=np.load(O/'h_cache.npy',mmap_mode='r'); U=np.load(O/'u_cache.npy',mmap_mode='r'); V=np.load(O/'v_cache.npy',mmap_mode='r')
valid_dem=np.loadtxt(P/'inputs_30m'/'dem.asc',skiprows=6)!=-9999
sec=pd.read_csv(P/'inputs_30m'/'analysis_sections.csv'); start=sec.iloc[0]; end=sec.iloc[-1]
releases=np.arange(10800,518400+1,10800,dtype=float)
if args.cohort_limit: releases=releases[:args.cohort_limit]
N=1000; B=len(releases)
px=np.full((B,N),np.nan); py=px.copy(); arrivals=px.copy()
state=np.zeros((B,N),dtype=np.int8) # 0 unreleased; 1 active; 2 arrived; 3 outside model; 4 timed out; 5 no source; 6 NoData
def spatial(a,x,y):
    return map_coordinates(a,[(top-y)/dx-.5,(x-left)/dx-.5],order=1,mode='constant',cval=0,prefilter=False)
def vel(t,x,y):
    f=min(int(t//600),len(U)-2); w=(t-f*600)/600
    return ((1-w)*spatial(U[f],x,y)+w*spatial(U[f+1],x,y),(1-w)*spatial(V[f],x,y)+w*spatial(V[f+1],x,y))
def signed(x,y): return (x-end.x)*end.normal_east+(y-end.y)*end.normal_north
next_release=0; t=0.; last_report=-1
ready_frame=-1
release_rows=[]
while t < releases[-1]+86400-1e-6:
    needed=min(int(t//600)+1,len(V)-1)
    if needed>ready_frame:
        # The cache writer writes h,u,v in frame order. Any written velocity in
        # the NEXT frame proves all fields in this needed frame were completed.
        deadline=time.monotonic()+7200
        while not (O/'cache_complete.txt').exists() and (needed+1>=len(V) or not np.any(V[needed+1])):
            if time.monotonic()>deadline: raise RuntimeError('Timed out waiting for velocity cache')
            time.sleep(2)
        ready_frame=needed
    while next_release<B and releases[next_release]<=t+1e-6:
        b=next_release; off=np.arange(-2000+3.75,2000,7.5)
        xx=start.x-off*start.normal_north; yy=start.y+off*start.normal_east
        uu,vv=vel(t,xx,yy); hh=spatial(H[int(round(t/600))],xx,yy)
        weights=np.maximum(uu*start.normal_east+vv*start.normal_north,0)*hh
        weights[hh<=.01]=0
        if weights.sum()<=0:
            state[b]=5
        else:
            cdf=np.cumsum(weights)/weights.sum()
            idx=np.searchsorted(cdf,(np.arange(N)+.5)/N)
            # Sub-bin stratification reduces identical trajectories from raster sections.
            lo=np.r_[0,cdf[:-1]][idx]; hi=cdf[idx]
            frac=((np.arange(N)+.5)/N-lo)/(hi-lo)
            loc=off[idx]+(frac-.5)*7.5
            px[b]=start.x-loc*start.normal_north+2*start.normal_east
            py[b]=start.y+loc*start.normal_east+2*start.normal_north
            state[b]=1
        release_rows.append(dict(release_s=t,source_positive_Q_approx_m3s=weights.sum()*7.5))
        next_release+=1
    expired=(state==1)&(t>=releases[:,None]+86400-1e-6)
    state[expired]=4
    b,j=np.where(state==1)
    dt=min(args.step,600-(t%600),releases[-1]+86400-t)
    if next_release<B: dt=min(dt,releases[next_release]-t)
    if len(b):
        xx=px[b,j]; yy=py[b,j]; u,v=vel(t,xx,yy)
        speed=np.hypot(u,v)
        dt=min(dt,dx*.4/max(float(speed.max()),.01))
        um,vm=vel(t+dt/2,xx+u*dt/2,yy+v*dt/2)
        xn=xx+um*dt; yn=yy+vm*dt
        old=signed(xx,yy); new=signed(xn,yn)
        lateral=(xn-end.x)*(-end.normal_north)+(yn-end.y)*end.normal_east
        hit=(old<0)&(new>=0)&(np.abs(lateral)<=2000)
        arrivals[b[hit],j[hit]]=t+dt*(-old[hit])/(new[hit]-old[hit])-releases[b[hit]]
        state[b[hit],j[hit]]=2
        out=(xn<left+dx/2)|(xn>left+(g['nx']-.5)*dx)|(yn<g['bottom']+dx/2)|(yn>top-dx/2)|~np.isfinite(xn+yn)
        state[b[out&~hit],j[out&~hit]]=3
        rc=np.clip(((top-yn)/dx).astype(int),0,g['ny']-1); cc=np.clip(((xn-left)/dx).astype(int),0,g['nx']-1)
        missing=~valid_dem[rc,cc]&~out&~hit
        state[b[missing],j[missing]]=6
        px[b,j]=xn; py[b,j]=yn
    t+=dt
    day=int(t//86400)
    if day!=last_report:
        print('particle simulation day',day,'active',int((state==1).sum()),'arrived',int((state==2).sum()),flush=True); last_report=day
state[state==1]=4
rows=[]
for b,r in enumerate(releases):
    a=arrivals[b]; finite=np.sort(a[np.isfinite(a)])
    row=dict(release_time=pd.Timestamp('2021-07-17')+pd.Timedelta(seconds=r),release_s=r,n_released=0 if np.all(state[b]==5) else N,T50_h=finite[499]/3600 if len(finite)>=500 else np.nan,T50_right_censored_24h=len(finite)<500,n_arrived_24h=len(finite),numerical_domain_exit_fraction=float(np.mean(state[b]==3)),numerical_nodata_fraction=float(np.mean(state[b]==6)),unarrived_in_domain_fraction=float(np.mean(state[b]==4)))
    for hours in [6,12,24]: row[f'F{hours}']=float(np.mean(a<=hours*3600)) if row['n_released'] else np.nan
    rows.append(row)
df=pd.DataFrame(rows)
inflow=pd.read_csv(P/'inputs_30m'/'inflow_hourly.csv'); inflow['time']=pd.to_datetime(inflow.time)
df=df.merge(inflow[['time','Q','HCI_fixed']],left_on='release_time',right_on='time',how='left').drop(columns='time')
suffix='' if args.step==10 and args.cohort_limit is None else f'_step{args.step:g}_n{B}'
df.to_csv(O/f'particle_cohorts{suffix}.csv',index=False)
pd.DataFrame(release_rows).to_csv(O/f'particle_release_flux{suffix}.csv',index=False)
np.savez_compressed(O/f'particle_outcomes{suffix}.npz',arrival_seconds=arrivals,state=state,release_seconds=releases,final_x=px,final_y=py)
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
axs[0,0].plot(df.release_time,df.T50_h,'o-',ms=3); axs[0,0].set_ylabel('T50 (h; all released particles)')
for k in ['F6','F12','F24']: axs[0,1].plot(df.release_time,df[k],label=k)
axs[0,1].set(ylabel='Arrival fraction',ylim=(-.03,1.03)); axs[0,1].legend()
sc=axs[1,0].scatter(df.HCI_fixed,df.T50_h,c=df.Q,cmap='viridis',s=25); fig.colorbar(sc,ax=axs[1,0],label='Q (m3/s)'); axs[1,0].set(xlabel='Fixed-RSE HCI',ylabel='T50 (h)')
axs[1,1].plot(df.release_time,df.numerical_domain_exit_fraction,label='Outside model domain'); axs[1,1].plot(df.release_time,df.unarrived_in_domain_fraction,label='In domain, not arrived by 24 h'); axs[1,1].set(ylabel='Fraction of all released',ylim=(-.03,1.03)); axs[1,1].legend(fontsize=8)
for ax in [axs[0,0],axs[0,1],axs[1,1]]: ax.tick_params(axis='x',rotation=30)
fig.suptitle('Exploratory passive advection; unvalidated bed and boundary assumptions')
for ext in ['png','pdf','svg']: fig.savefig(O/f'particle_diagnostics{suffix}.{ext}',dpi=300)
print(df.to_string(index=False),flush=True)
if default_run: (O/'particle_complete.txt').write_text('All default cohorts and figures complete.\n')
