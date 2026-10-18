"""Generate compact key HCI and validation figures from authorized source tables."""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('--hci',type=Path,required=True)
parser.add_argument('--matched-cohorts',type=Path,required=True)
parser.add_argument('--output',type=Path,default=ROOT/'figures')
args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
colors={"HYK":"#1f77b4","JHT":"#e66101","GC":"#2ca02c"}
h=pd.read_csv(args.hci,parse_dates=["time"])
fig,ax=plt.subplots(figsize=(8,3.8))
for code,g in h.groupby("station_code"):
    ax.plot(g.time,g.HCI_ref,label=code,color=colors[code],lw=1.7)
ax.set(ylabel="HCI",xlabel="Time (2021)",ylim=(0,1));ax.legend(frameon=False,ncol=3);fig.tight_layout();fig.savefig(args.output/"reproduced_HCI_process.png",dpi=300);plt.close(fig)
m=pd.read_csv(args.matched_cohorts)
fig,ax=plt.subplots(figsize=(5,4))
for code,g in m[m.scenario.isin(colors)].groupby("scenario"):
    ax.scatter(g.HCI_ref,g.T50_h,s=18,label=code,color=colors[code])
ax.set(xlabel="HCI at release",ylabel="$T_{50}$ (h)");ax.legend(frameon=False);fig.tight_layout();fig.savefig(args.output/"reproduced_HCI_T50.png",dpi=300);plt.close(fig)
print(args.output)
