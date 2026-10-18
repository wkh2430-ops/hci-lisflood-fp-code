"""Recalculate Cq(t) and HCI(t) without changing the published scientific logic."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True,help="CSV containing time, station_code, and Q_m3_s or Q_ref")
    ap.add_argument("--config",default="configs/stations.json")
    ap.add_argument("--output",default="results/HCI_hourly_recalculated.csv")
    a=ap.parse_args(); root=Path(__file__).resolve().parents[2]
    cfg=json.loads((root/a.config).read_text(encoding="utf-8")); d=pd.read_csv(root/a.input)
    d["time"]=pd.to_datetime(d["time"])
    qcol="Q_ref" if "Q_ref" in d else "Q_m3_s"
    required={"time","station_code",qcol}
    missing=required-set(d.columns)
    if missing: raise ValueError(f"Missing required columns: {sorted(missing)}")
    rows=[]
    for code,g in d.groupby("station_code",sort=False):
        s=cfg["stations"][code]; x=g.copy(); q=x[qcol].astype(float)
        x["Qbank_m3_s"]=s["Qbank_m3_s"]
        # The final 2021-07-17 output used rounded station-specific k in
        # Cq=1/{1+exp[-k(Q-Qbank)]}. This is retained for exact reproduction.
        x["Cq_recalculated"]=1/(1+np.exp(-s["legacy_k_s_m3"]*(q-s["Qbank_m3_s"])))
        x["RSE_fixed_recalculated"]=s["RSE_fixed"]
        x["HCI_recalculated"]=x.Cq_recalculated*s["RSE_fixed"]
        rows.append(x)
    out=pd.concat(rows).sort_values(["station_code","time"])
    p=root/a.output;p.parent.mkdir(parents=True,exist_ok=True);out.to_csv(p,index=False)
    if "Cq_ref" in out:
        print("max_abs_Cq_difference",float(np.max(np.abs(out.Cq_recalculated-out.Cq_ref))))
    if "HCI_ref" in out:
        print("max_abs_HCI_difference",float(np.max(np.abs(out.HCI_recalculated-out.HCI_ref))))
    print(p)
if __name__=="__main__":main()
