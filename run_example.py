"""Run and verify the public synthetic HCI example."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
cfg = json.loads((ROOT / "configs" / "stations.json").read_text(encoding="utf-8"))
data = pd.read_csv(ROOT / "examples" / "synthetic_hci_input_and_expected.csv")

max_cq_difference = 0.0
max_hci_difference = 0.0
for code, group in data.groupby("station_code"):
    station = cfg["stations"][code]
    q = group["Q_m3_s"].to_numpy(float)
    cq = 1 / (1 + np.exp(-station["legacy_k_s_m3"] * (q - station["Qbank_m3_s"])))
    hci = station["RSE_fixed"] * cq
    max_cq_difference = max(max_cq_difference, float(np.max(np.abs(cq - group["Cq"]))))
    max_hci_difference = max(max_hci_difference, float(np.max(np.abs(hci - group["HCI"]))))
    np.testing.assert_allclose(cq, group["Cq"], rtol=0, atol=1e-12)
    np.testing.assert_allclose(hci, group["HCI"], rtol=0, atol=1e-12)

print(f"Synthetic HCI example reproduced (max |ΔCq|={max_cq_difference:.3g}, max |ΔHCI|={max_hci_difference:.3g}).")
