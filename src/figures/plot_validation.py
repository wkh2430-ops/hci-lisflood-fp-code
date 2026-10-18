"""Regenerate the manuscript validation figures from prepared source tables."""
from pathlib import Path
import argparse
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.patches import Patch
parser=argparse.ArgumentParser()
parser.add_argument('--source-data',type=Path,default=Path('data/processed/validation'))
parser.add_argument('--stage-data',type=Path,default=Path('data/processed/hci/flood_stage_boundaries.csv'))
parser.add_argument('--output',type=Path,default=Path('figures'))
args=parser.parse_args()
D=args.source_data;F=args.output;F.mkdir(parents=True,exist_ok=True)
x=pd.read_csv(D/'matched_cohorts.csv',parse_dates=['time']);h=pd.read_csv(D/'HCI_hourly.csv',parse_dates=['time']);s=pd.read_csv(D/'summary.csv').set_index('scenario');cdf=pd.read_csv(D/'cdf.csv');rep=pd.read_csv(D/'representative_cohorts.csv');codes=['HYK','JHT','GC'];colors=dict(zip(codes,['#1f77b4','#ff7f0e','#2ca02c']))
stage=pd.read_csv(args.stage_data).set_index('station_code');shades=['#e8f3e8','#fff6df','#fff0df','#e8eef8'];labels=['Low-flow stage','Rising stage','Peak stage','Recession stage'];markers=['o','s','^']
plt.rcParams.update({'font.family':['Times New Roman'],'mathtext.fontset':'stix','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':1,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.unicode_minus':False})
def save(fig,name):
 import re
 from matplotlib.text import Text
 for ax in fig.axes:
  for loc in ('left','center','right'):
   title=ax.get_title(loc=loc)
   if not title:
    continue
   panel=re.match(r'^(\([a-z]\))',title)
   station=re.search(r'\b(HYK|JHT|GC)\b',title)
   if name.startswith(('Fig_V1_','Fig_V2_')) and panel and station:
    rho=re.search(r'ρ = -?[0-9.]+',title)
    clean=panel.group(1)+' '+station.group(1)
    if rho:
     clean+=': '+rho.group(0)
   else:
    clean=panel.group(1) if panel else ''
   ax.set_title(clean,loc=loc)
 text_items=[t.get_text() for t in fig.findobj(Text)]
 assert not any(any(n in t for n in ('Huayuankou','Jiahetan','Gaocun')) for t in text_items)
 assert not any(re.search(r'[\u3400-\u9fff]',t) for t in text_items), text_items
 assert not any('Gaocun (630' in t or '高村' in t for t in text_items)

 for ext in ['png','pdf','svg','tiff']:
  kw={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
  fig.savefig(F/(name+'.'+ext),dpi=300,bbox_inches='tight',**kw)
 plt.close(fig)
def grid(a):a.grid(axis='y',ls='--',alpha=.22);a.set_axisbelow(True)
def date(a):a.xaxis.set_major_locator(md.DayLocator());a.xaxis.set_major_formatter(md.DateFormatter('%m-%d'));a.set_xlim(pd.Timestamp('2021-07-17'),pd.Timestamp('2021-07-24'))
def shade(a,code):
 r=stage.loc[code];bounds=[pd.Timestamp('2021-07-17')]+[pd.Timestamp(r[k]) for k in ['rising_start','peak_stage_start','recession_start']]+[pd.Timestamp('2021-07-24')]
 for l,r,c in zip(bounds[:-1],bounds[1:],shades):a.axvspan(l,r,color=c,alpha=.5,lw=0)
# V1: exact release-time pairing, not shifted downstream arrival peaks.
fig,axs=plt.subplots(3,1,figsize=(8,8),sharex=True,layout='constrained')
for i,(a,code) in enumerate(zip(axs,codes)):
 d=x[x.scenario==code];r=h[h.station_code==code];shade(a,code);grid(a);a.plot(r.time,r.HCI_ref,color=colors[code],lw=1.8,label='HCI(t)');a.axhline(s.loc[code,'RSE'],c='gray',ls='--',lw=1,label='RSE');a.set_ylim(0,1.02);a.set_ylabel('HCI / RSE');b=a.twinx();b.spines['right'].set_visible(True);b.plot(d.time,d.T50_h,'o-',c='#333333',ms=3,lw=1.2,label='T₅₀');b.set_ylim(1.5,4.8);b.set_ylabel('T₅₀ (h)');a.set_title(f'({chr(97+i)}) {code}',loc='left');date(a)
 if i==0:
  ln=a.lines+b.lines;a.legend(ln,[l.get_label() for l in ln],loc='upper left',ncol=3,fontsize=9)
axs[-1].set_xlabel('Particle release time / HCI time (2021)');fig.legend(handles=[Patch(color=c,alpha=.5,label=l) for c,l in zip(shades,labels)],loc='outside upper center',ncol=4,fontsize=9);save(fig,'Fig_V1_HCI_T50_Time_Series')
# V2 relationships, no pooled regression or p values.
fig,axs=plt.subplots(2,3,figsize=(10,6),layout='constrained')
for j,code in enumerate(codes):
 d=x[x.scenario==code]
 for i,(metric,ylab) in enumerate([('T50_h','T₅₀ (h)'),('F24','F₂₄ (%)')]):
  a=axs[i,j]
  for phase,mk in [('Pre-peak','o'),('Post-peak','^')]:
   dd=d[d.phase==phase];a.scatter(dd.HCI_ref,dd[metric]*(100 if i else 1),s=22,marker=mk,c=colors[code],alpha=.8,label=phase)
  grid(a);a.set_xlabel('HCI');a.set_ylabel(ylab);a.set_ylim((89,100) if i else (1.5,4.8));rho=s.loc[code,'rho_HCI_F24' if i else 'rho_HCI_T50'];a.set_title(f'({chr(97+i*3+j)}) {code}: ρ = {rho:.3f}');a.legend(fontsize=8)
save(fig,'Fig_V2_HCI_Transport_Relationships')
# V3 unconditional cumulative arrival, all released particles denominator.
fig,axs=plt.subplots(2,2,figsize=(9,6.6),layout='constrained')
for i,(a,phase,title) in enumerate(zip(axs.flat,['All cohorts','Pre-peak 24 h','Near peak','Post-peak 24 h'],['All 48 cohorts pooled equally','Released ~24 h before local peak','Released near local peak','Released ~24 h after local peak'])):
 for code,mk in zip(codes,markers):
  d=cdf[(cdf.scenario==code)&(cdf.stage==phase)];a.plot(d.elapsed_h,d.F*100,color=colors[code],lw=1.7,label=code,marker=mk,markevery=80,ms=3)
 a.axhline(50,c='gray',ls=':',lw=.8);a.axhline(90,c='gray',ls=':',lw=.8);a.set(xlim=(0,24),ylim=(0,101),xlabel='Time since release (h)',ylabel='Cumulative arrival (%)');a.set_title(f'({chr(97+i)}) {title}',loc='left');grid(a);a.legend(fontsize=9)
save(fig,'Fig_V3_Cumulative_Arrival')
# V4 summary, distinguishes median of cohorts vs pooled quantiles.
fig,axs=plt.subplots(2,2,figsize=(9,6.5),layout='constrained');pos=np.arange(3);cc=[colors[c] for c in codes];main=s.loc[codes]
a=axs[0,0];a.bar(pos-.18,main.HCI_mean_hourly,.34,color=cc,label='Full-event hourly mean');a.bar(pos+.18,main.HCI_mean_releases,.34,color=cc,alpha=.35,hatch='//',label='Mean at 48 release times');a.set_ylabel('Mean HCI');a.set_ylim(0,.48);a.legend(fontsize=8);a.set_title('(a) HCI means across sampling windows',loc='left')
a=axs[0,1]
for j,c in enumerate(codes):
 vals=x[x.scenario==c].T50_h;lo,med,hi=np.percentile(vals,[25,50,75]);a.errorbar(j,med,yerr=[[med-lo],[hi-med]],fmt=markers[j],c=colors[c],capsize=5,ms=7);a.text(j,hi+.10,f'{med:.2f}',ha='center')
a.set_ylabel('Cohort T₅₀: median and IQR (h)');a.set_ylim(1.5,5);a.set_title('(b) Median cohort travel time',loc='left')
a=axs[1,0]
for j,c in enumerate(codes):
 a.plot([j-.14,j+.14],[s.loc[c,'T50_pooled'],s.loc[c,'T90_pooled']],c=colors[c]);a.scatter(j-.14,s.loc[c,'T50_pooled'],marker='o',c=colors[c]);a.scatter(j+.14,s.loc[c,'T90_pooled'],marker='s',c=colors[c])
a.plot([],[],'ko',label='Pooled T₅₀');a.plot([],[],'ks',label='Pooled T₉₀');a.legend(fontsize=8);a.set_ylabel('Pooled arrival-time quantile (h)');a.set_ylim(0,6.5);a.set_title('(c) Pooled travel-time quantiles',loc='left')
a=axs[1,1]
for j,c in enumerate(codes):
 d=x[x.scenario==c];mean=d.F24.mean()*100;a.errorbar(j,mean,yerr=[[mean-d.F24.min()*100],[d.F24.max()*100-mean]],fmt=markers[j],c=colors[c],capsize=5,ms=7);a.text(j,mean+.35,f'{mean:.2f}%',ha='center',fontsize=9)
a.set_ylim(88,101);a.set_ylabel('F₂₄: mean and cohort range (%)');a.set_title('(d) 24 h arrival fraction',loc='left')
for a in axs.flat:a.set_xticks(pos,codes);grid(a)
save(fig,'Fig_V4_Three_Reach_Comparison')
# V5 common flow dependence, paired descriptive coefficients without p values.
fig,axs=plt.subplots(1,2,figsize=(9,3.6),layout='constrained');diagnostics=[]
for j,c in enumerate(codes):
 d=x[x.scenario==c]
 for i,diff in enumerate([False,True]):
  vals=d[['HCI_ref','Q_ref','T50_h']].diff().dropna() if diff else d[['HCI_ref','Q_ref','T50_h']]
  rh=vals.HCI_ref.corr(vals.T50_h,method='spearman');rq=vals.Q_ref.corr(vals.T50_h,method='spearman');diagnostics.append(dict(station=c,first_difference=diff,rho_HCI_T50=rh,rho_Q_T50=rq))
  axs[i].bar(j-.16,rh,.30,color=colors[c]);axs[i].bar(j+.16,rq,.30,color=colors[c],alpha=.35,hatch='//')
for i,a in enumerate(axs):
 a.set(xticks=pos,xticklabels=codes,ylim=(-1.05,0),ylabel='Spearman ρ');a.set_title(['(a) Original series: HCI and Q share ranks','(b) First differences: trend diagnostic'][i],loc='left',fontsize=10);grid(a);a.legend(handles=[Patch(facecolor='gray',label='HCI–T₅₀'),Patch(facecolor='gray',alpha=.35,hatch='//',label='Q–T₅₀')],loc='upper right',fontsize=8)
save(fig,'Fig_V5_Shared_Discharge_Diagnostic')
print('Five figure groups exported')
