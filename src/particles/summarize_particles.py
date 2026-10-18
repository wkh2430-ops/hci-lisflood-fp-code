"""Recompute cohort T50/T90/F6/F12/F24 using all released particles as denominator."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def quantile_all_released(v,q):
    vals=np.sort(v[np.isfinite(v)]); needed=int(np.ceil(q*v.size))
    return vals[needed-1]/3600 if len(vals)>=needed else np.nan
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--input",type=Path,required=True,help="particle_outcomes.npz");ap.add_argument("--scenario",required=True);ap.add_argument("--output",type=Path);a=ap.parse_args()
    root=Path(__file__).resolve().parents[2]
    z=np.load(a.input);arr=z["arrival_seconds"]
    rows=[]
    for b,v in enumerate(arr):
        row={"scenario":a.scenario,"batch":b,"release_s":float(z["release_seconds"][b]),"n_released":int(v.size),"n_arrived_24h":int(np.isfinite(v).sum()),"T50_h":quantile_all_released(v,.5),"T90_h":quantile_all_released(v,.9)}
        for h in [6,12,24]:row[f"F{h}"]=float(np.mean(v<=h*3600))
        rows.append(row)
    out=pd.DataFrame(rows);p=a.output or root/"results"/f"{a.scenario}_particle_metrics_recalculated.csv";p.parent.mkdir(parents=True,exist_ok=True);out.to_csv(p,index=False);print(p)
if __name__=="__main__":main()
