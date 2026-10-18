"""Fixed-RSE sensitivity analysis for Qbank, beta, and activation interval."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def logistic(q: pd.Series, qbank: float, beta: float) -> np.ndarray:
    z = np.clip(beta * (q.to_numpy(dtype=float) / qbank - 1.0), -700, 700)
    return 1.0 / (1.0 + np.exp(-z))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hci", type=Path, default=Path("data/processed/hci/HCI_hourly.csv"))
    parser.add_argument("--config", type=Path, default=Path("configs/stations.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/sensitivity"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    hci = pd.read_csv(args.hci, parse_dates=["time"])
    config = json.loads(args.config.read_text(encoding="utf-8"))
    beta0 = float(config["formula"]["beta"])
    qbank_factors = config["sensitivity"]["Qbank_factors"]
    beta_values = config["sensitivity"]["beta_values"]
    ranges = config["sensitivity"]["activation_ranges_Qbank"]

    rows = []
    for code, group in hci.groupby("station_code", sort=True):
        station = config["stations"][code]
        qbank0 = float(station["Qbank_m3_s"])
        rse = float(station["RSE_fixed"])
        for factor in qbank_factors:
            cq = logistic(group.Q_m3_s, qbank0 * factor, beta0)
            for time, q, c in zip(group.time, group.Q_m3_s, cq):
                rows.append((time, code, "Qbank", f"{factor:.1f}Qbank", q, qbank0 * factor, beta0, rse, c, rse * c))
        for beta in beta_values:
            cq = logistic(group.Q_m3_s, qbank0, float(beta))
            for time, q, c in zip(group.time, group.Q_m3_s, cq):
                rows.append((time, code, "beta", f"beta={beta:.8g}", q, qbank0, float(beta), rse, c, rse * c))
        for low, high in ranges:
            beta = float(np.log(9.0) / (1.0 - low))
            cq = logistic(group.Q_m3_s, qbank0, beta)
            label = f"Cq0.1-0.9={low:.2f}-{high:.2f}Qbank"
            for time, q, c in zip(group.time, group.Q_m3_s, cq):
                rows.append((time, code, "activation_range", label, q, qbank0, beta, rse, c, rse * c))

    columns = ["time", "station_code", "parameter", "scenario", "Q_m3_s", "Qbank_scenario_m3_s", "beta", "RSE_fixed", "Cq", "HCI"]
    hourly = pd.DataFrame(rows, columns=columns)
    hourly.to_csv(args.output / "HCI_sensitivity_hourly.csv", index=False)
    baseline_means = hci.groupby("station_code").HCI.mean()
    summary = (
        hourly.groupby(["station_code", "parameter", "scenario"], sort=False)
        .agg(HCI_mean=("HCI", "mean"), HCI_max=("HCI", "max"), Cq_mean=("Cq", "mean"), Cq_max=("Cq", "max"))
        .reset_index()
    )
    summary["HCI_mean_change_pct_vs_exact_baseline"] = summary.apply(
        lambda row: 100 * (row.HCI_mean / baseline_means[row.station_code] - 1), axis=1
    )
    summary.to_csv(args.output / "HCI_sensitivity_summary.csv", index=False)
    print(f"Generated fixed-RSE sensitivity tables in {args.output}")


if __name__ == "__main__":
    main()
