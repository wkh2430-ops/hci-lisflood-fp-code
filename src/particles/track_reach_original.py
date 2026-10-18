"""Passive particle diagnostic with explicit wet/dry and numerical exit states.
Run only after event water-balance checks; no HCI-dependent hydraulic parameters.
"""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from scipy.ndimage import map_coordinates
ap=argparse.ArgumentParser();ap.add_argument('reach');ap.add_argument('--step',type=float,default=10);ap.add_argument('--limit',type=int);a=ap.parse_args()
P=Path(a.reach);O=P/'analysis';O.mkdir(exist_ok=True);R=P/'event/results'
g=json.loads((P/'manifest.json').read_text())['grid'];dx=g['dx'];left=g['left'];top=g['top'];bottom=g['bottom'];nx=g['nx'];ny=g['ny']
z=np.loadtxt(P/'inputs/dem.asc',skiprows=6);valid=np.isfinite(z)&(z!=-9999)
m=pd.read_csv(R/'run.mass',sep=r'\s+');files=sorted(R.glob('run-*.wd'))
if len(files)!=1009 or m.Time.iloc[-1]<604800:raise RuntimeError('Incomplete event')
times=m.Time.to_numpy()[np.searchsorted(m.Time.to_numpy(),np.arange(1009)*600)]
if np.any(np.diff(times)<=0) or np.max(times-np.arange(1009)*600)>m.Tstep.max()+.02:raise RuntimeError('Frame-time recovery failed')
pd.DataFrame({'frame':range(1009),'actual_seconds':times}).to_csv(O/'frame_times.csv',index=False)
if not (O/'cache_complete.txt').exists():
 caches=[np.lib.format.open_memmap(O/(k+'_cache.npy'),mode='w+',dtype='float32',shape=(1009,ny,nx)) for k in ['h','u','v']]
 for i,f in enumerate(files):
  h=np.loadtxt(f,skiprows=6,dtype='float32');vx=np.loadtxt(f.with_suffix('.Vx'),skiprows=6,dtype='float32');vy=np.loadtxt(f.with_suffix('.Vy'),skiprows=6,dtype='float32')
  vx[vx<-999]=0;vy[vy<-999]=0;u=(vx[:,:-1]+vx[:,1:])/2;v=-(vy[:-1]+vy[1:])/2
  for cache,field in zip(caches,[h,u,v]):
   field[~valid]=0
   if not np.isfinite(field).all():raise RuntimeError('Nonfinite field')
   cache[i]=field
  if i%100==0:print('cache',i,flush=True)
 for ca in caches:ca.flush()
 (O/'cache_complete.txt').write_text('1009 checked frames')
H,U,V=[np.load(O/(k+'_cache.npy'),mmap_mode='r') for k in ['h','u','v']]
sec=pd.read_csv(P/'inputs/sections_local.csv');start=sec.iloc[0];end=sec.iloc[-1]
rels=np.arange(10800,518401,10800,dtype=float)
if a.limit:rels=rels[:a.limit]
B=len(rels);N=1000;px=np.full((B,N),np.nan);py=px.copy();arrival=px.copy();state=np.zeros((B,N),'int8');drytime=np.zeros((B,N));blocked=np.zeros((B,N));trajectory=[]
def spatial(f,x,y,order=1):return map_coordinates(f,[(top-y)/dx-.5,(x-left)/dx-.5],order=order,mode='constant',cval=0,prefilter=False)
def fields(t,x,y):
 i=int(np.clip(np.searchsorted(times,t,side='right')-1,0,1007));w=np.clip((t-times[i])/(times[i+1]-times[i]),0,1)
 return tuple((1-w)*spatial(f[i],x,y)+w*spatial(f[i+1],x,y) for f in [H,U,V])
def inside(x,y):return (x>=left+dx/2)&(x<=left+(nx-.5)*dx)&(y>=bottom+dx/2)&(y<=top-dx/2)&np.isfinite(x+y)
def signed(x,y):return (x-end.x)*end.normal_east+(y-end.y)*end.normal_north
t=rels[0];nextbatch=0;lastframe=-1;release_rows=[]
while t<rels[-1]+86400-1e-7:
 if nextbatch<B and abs(t-rels[nextbatch])<1e-5:
  b=nextbatch;off=np.arange(-1996.25,2000,7.5);xx=start.x-off*start.normal_north;yy=start.y+off*start.normal_east;hh,uu,vv=fields(t,xx,yy)
  weights=np.maximum(uu*start.normal_east+vv*start.normal_north,0)*hh
  fullbin=inside(xx-3.75*start.normal_north+2*start.normal_east,yy+3.75*start.normal_east+2*start.normal_north)&inside(xx+3.75*start.normal_north+2*start.normal_east,yy-3.75*start.normal_east+2*start.normal_north)
  weights[(hh<.01)|~fullbin]=0
  if weights.sum()<=0:state[b]=5
  else:
   cdf=np.cumsum(weights)/weights.sum();q=(np.arange(N)+.5)/N;ix=np.searchsorted(cdf,q);lo=np.r_[0,cdf[:-1]][ix];frac=(q-lo)/(cdf[ix]-lo);loc=off[ix]+(frac-.5)*7.5
   px[b]=start.x-loc*start.normal_north+2*start.normal_east;py[b]=start.y+loc*start.normal_east+2*start.normal_north;state[b]=1
  release_rows.append({'release_s':t,'positive_Q_approx_m3s':weights.sum()*7.5});nextbatch+=1
 state[(state==1)&(t>=rels[:,None]+86400-1e-6)]=4
 dt=min(a.step,rels[-1]+86400-t)
 if nextbatch<B:dt=min(dt,rels[nextbatch]-t)
 active_end=rels+86400;future=active_end[active_end>t+1e-6]
 if len(future):dt=min(dt,future.min()-t)
 b,j=np.where(state==1)
 if len(b):
  xx=px[b,j];yy=py[b,j];h,u,v=fields(t,xx,yy);dry=h<.01;u[dry]=0;v[dry]=0;dt=min(dt,dx*.4/max(float(np.hypot(u,v).max()),.01))
  hm,um,vm=fields(t+dt/2,xx+u*dt/2,yy+v*dt/2);um[hm<.01]=0;vm[hm<.01]=0
  xn=xx+dt*um;yn=yy+dt*vm;out=~inside(xn,yn);missing=(spatial(valid.astype('uint8'),xn,yn,0)<.5)&~out;hn,_,_=fields(t+dt,xn,yn);barrier=((hn<.01)|(hm<.01))&~out&~missing;xn[barrier]=xx[barrier];yn[barrier]=yy[barrier]
  drytime[b,j]+=dt*dry;blocked[b,j]+=dt*barrier
  old=signed(xx,yy);new=signed(xn,yn);lateral=-(xn-end.x)*end.normal_north+(yn-end.y)*end.normal_east
  hit=(old<0)&(new>=0)&(abs(lateral)<=2000)&~out&~missing
  arrival[b[hit],j[hit]]=t+dt*(-old[hit])/(new[hit]-old[hit])-rels[b[hit]];state[b[hit],j[hit]]=2;state[b[out],j[out]]=3;state[b[missing],j[missing]]=6;px[b,j]=xn;py[b,j]=yn
 t+=dt;fr=int(t//3600)
 if fr!=lastframe:
  for b in range(B):
   if t>=rels[b]:
    for j in range(0,N,50):trajectory.append((t,b,j,px[b,j],py[b,j],state[b,j]))
  print('hour',fr,'arrived',int((state==2).sum()),flush=True);lastframe=fr
state[state==1]=4;rows=[]
for b,r in enumerate(rels):
 v=np.sort(arrival[b,np.isfinite(arrival[b])]);row=dict(release_time=pd.Timestamp('2021-07-17')+pd.Timedelta(seconds=r),release_s=r,n_released=0 if np.all(state[b]==5) else N,T50_h=v[499]/3600 if len(v)>=500 else np.nan,T50_right_censored_24h=len(v)<500,n_arrived_24h=len(v),numerical_domain_exit_fraction=float(np.mean(state[b]==3)),numerical_nodata_fraction=float(np.mean(state[b]==6)),unarrived_in_domain_fraction=float(np.mean(state[b]==4)),any_dry_exposure_fraction=float(np.mean(drytime[b]>0)),any_dry_boundary_block_fraction=float(np.mean(blocked[b]>0)))
 for hours in [6,12,24]:row[f'F{hours}']=float(np.mean(arrival[b]<=hours*3600)) if row['n_released'] else np.nan
 rows.append(row)
df=pd.DataFrame(rows);hyd=pd.read_csv(P/'inputs/inflow_hourly.csv');hyd.time=pd.to_datetime(hyd.time);df=df.merge(hyd,left_on='release_time',right_on='time',validate='one_to_one')
suffix='' if not a.limit and a.step==10 else f'_step{a.step:g}_n{B}'
df.to_csv(O/f'particle_cohorts{suffix}.csv',index=False,encoding='utf-8-sig');np.savez_compressed(O/f'particle_outcomes{suffix}.npz',arrival_seconds=arrival,state=state,release_seconds=rels,final_x=px,final_y=py,dry_seconds=drytime,dry_boundary_block_seconds=blocked)
pd.DataFrame(trajectory,columns=['seconds','batch','particle','x','y','state']).to_csv(O/f'trajectory_samples{suffix}.csv',index=False);pd.DataFrame(release_rows).to_csv(O/f'release_flux{suffix}.csv',index=False)
print('COMPLETE',str(O),flush=True)
