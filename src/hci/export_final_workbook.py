"""Export normalized source tables from the paper-final fixed-RSE workbook."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("data/processed/hci"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    hourly = pd.read_excel(args.workbook, sheet_name=0)
    fixed = pd.read_excel(args.workbook, sheet_name=1)
    cq = pd.read_excel(args.workbook, sheet_name=2)
    station = pd.read_excel(args.workbook, sheet_name=3)
    stages = pd.read_excel(args.workbook, sheet_name=4)

    out = pd.DataFrame(
        {
            "time": pd.to_datetime(hourly.iloc[:, 0]),
            "station_code": hourly["station_code"].astype(str).str.strip(),
            "Q_m3_s": pd.to_numeric(hourly["Q(t)(m3/s)"], errors="raise"),
            "RSE_fixed": pd.to_numeric(hourly["RSE_fixed_20210717_0000"], errors="raise"),
            "Cq": pd.to_numeric(hourly["Cq(t)=1/(1+exp{-k[Q(t)-Qc]})"], errors="raise"),
            "HCI": pd.to_numeric(hourly["HCI_0717RSE=RSE_fixed*Cq(t)"], errors="raise"),
        }
    ).sort_values(["station_code", "time"])
    if len(out) != 507 or out.groupby("station_code").size().to_dict() != {"GC": 169, "HYK": 169, "JHT": 169}:
        raise ValueError("Unexpected final-workbook record count")
    out.to_csv(args.output / "HCI_hourly.csv", index=False, encoding="utf-8")

    fixed_out = pd.DataFrame(
        {
            "station_code": fixed["station_code"].astype(str).str.strip(),
            "reference_time": pd.to_datetime(fixed.iloc[:, 2]),
            "Q_reference_m3_s": pd.to_numeric(fixed.iloc[:, 3]),
            "water_level_reference_m": pd.to_numeric(fixed.iloc[:, 4]),
            "water_depth_reference_m": pd.to_numeric(fixed.iloc[:, 5]),
            "RSE_fixed": pd.to_numeric(fixed["RSE_0717_0000"]),
            "relative_delay": pd.to_numeric(fixed["Dr=1/RSE-1"]),
        }
    )
    fixed_out.to_csv(args.output / "fixed_RSE_reference.csv", index=False, encoding="utf-8")

    cq_out = pd.DataFrame(
        {
            "station_code": cq["station_code"].astype(str).str.strip(),
            "Qbank_m3_s": pd.to_numeric(cq["Qbank_Qc(m3/s)"]),
            "Ql_0.75Qbank_m3_s": pd.to_numeric(cq["Q10=0.75Qbank(m3/s)"]),
            "Qh_1.25Qbank_m3_s": pd.to_numeric(cq["Q90=1.25Qbank(m3/s)"]),
            "legacy_k_s_m3": pd.to_numeric(cq["k"]),
        }
    )
    cq_out.to_csv(args.output / "Cq_parameters.csv", index=False, encoding="utf-8")

    station_out = pd.DataFrame(
        {
            "station_code": station["station_code"].astype(str).str.strip(),
            "n_hourly_records": pd.to_numeric(station.iloc[:, 2]),
            "reference_time": pd.to_datetime(station.iloc[:, 3]),
            "RSE_fixed": pd.to_numeric(station.iloc[:, 4]),
            "relative_delay": pd.to_numeric(station.iloc[:, 5]),
            "Cq_min": pd.to_numeric(station.iloc[:, 6]),
            "Cq_max": pd.to_numeric(station.iloc[:, 7]),
            "Cq_mean": pd.to_numeric(station.iloc[:, 8]),
            "HCI_min": pd.to_numeric(station.iloc[:, 9]),
            "HCI_max": pd.to_numeric(station.iloc[:, 10]),
            "HCI_mean": pd.to_numeric(station.iloc[:, 11]),
            "HCI_peak_time": pd.to_datetime(station.iloc[:, 12]),
            "Q_at_HCI_peak_m3_s": pd.to_numeric(station.iloc[:, 13]),
            "Cq_at_HCI_peak": pd.to_numeric(station.iloc[:, 14]),
        }
    )
    station_out.to_csv(args.output / "HCI_station_summary.csv", index=False, encoding="utf-8")

    stage_out = pd.DataFrame(
        {
            "station_code": stages["station_code"].astype(str).str.strip(),
            "Qbase_m3_s": pd.to_numeric(stages.iloc[:, 2]),
            "Qpeak_m3_s": pd.to_numeric(stages.iloc[:, 3]),
            "Q10_event_m3_s": pd.to_numeric(stages.iloc[:, 4]),
            "Q90_event_m3_s": pd.to_numeric(stages.iloc[:, 5]),
            "rising_start": pd.to_datetime(stages.iloc[:, 6]),
            "peak_stage_start": pd.to_datetime(stages.iloc[:, 7]),
            "peak_time": pd.to_datetime(stages.iloc[:, 8]),
            "recession_start": pd.to_datetime(stages.iloc[:, 9]),
        }
    )
    stage_out.to_csv(args.output / "flood_stage_boundaries.csv", index=False, encoding="utf-8")
    print(f"Exported final HCI tables to {args.output}")


if __name__ == "__main__":
    main()
