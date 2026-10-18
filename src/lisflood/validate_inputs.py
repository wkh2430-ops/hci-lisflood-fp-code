"""Validate packaged .par files and report required external inputs; does not run LISFLOOD-FP."""
from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parents[2]
keys={"DEMfile","manningfile","startfile","bcifile","bdyfile","sim_time","saveint","massint","initial_tstep","cfl"}
rows=[]
for p in sorted((ROOT/"configs/lisflood").glob("*.par")):
    vals={}
    for line in p.read_text(errors="replace").splitlines():
        part=line.split()
        if part and part[0] in keys:vals[part[0]]=" ".join(part[1:])
    rows.append({"file":str(p.relative_to(ROOT)),"parameters":vals,"configuration_complete":all(k in vals for k in ["DEMfile","manningfile","startfile","bcifile","bdyfile"]),"external_inputs_present_in_release":False})
out=ROOT/"results/lisflood_config_check.json";out.write_text(json.dumps(rows,indent=2),encoding="utf-8");print(out)
