"""Build the particle-validation source tables used by the manuscript figures."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SCENARIOS = ("HYK", "JHT", "GC")


def quantile_all_released(values: np.ndarray, q: float) -> float:
    finite = np.sort(values[np.isfinite(values)])
    required = int(np.ceil(q * values.size))
    return float(finite[required - 1] / 3600.0) if len(finite) >= required else np.nan


def duration_above(values: np.ndarray, threshold: float) -> float:
    total = 0.0
    for left, right in zip(values[:-1], values[1:]):
        if left >= threshold and right >= threshold:
            total += 1.0
        elif (left >= threshold) != (right >= threshold):
            total += (left - threshold) / (left - right) if left >= threshold else (right - threshold) / (right - left)
    return float(total)


def spearman(x: pd.Series, y: pd.Series) -> float:
    pair = pd.concat([x, y], axis=1).dropna()
    return float(pair.iloc[:, 0].rank(method="average").corr(pair.iloc[:, 1].rank(method="average")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hci", type=Path, default=Path("data/processed/hci/HCI_hourly.csv"))
    parser.add_argument("--particles", type=Path, default=Path("data/processed/particles"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/validation"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    hci = pd.read_csv(args.hci, parse_dates=["time"]).rename(
        columns={"Q_m3_s": "Q_ref", "HCI": "HCI_ref", "RSE_fixed": "RSE_ref", "Cq": "Cq_ref"}
    )
    hci.to_csv(args.output / "HCI_hourly.csv", index=False)
    summaries, matched, cdf, representative, provenance = [], [], [], [], []

    for code in SCENARIOS:
        station = "GC" if code.startswith("GC") else code
        reference = hci[hci.station_code == station].sort_values("time")
        cohort_path = args.particles / code / "particle_cohorts.csv"
        outcome_path = args.particles / code / "particle_outcomes.npz"
        cohort = pd.read_csv(cohort_path)
        cohort["time"] = pd.to_datetime(cohort["release_time"])
        cohort = cohort.merge(reference[["time", "Q_ref", "HCI_ref", "RSE_ref", "Cq_ref"]], on="time", validate="one_to_one")
        arrays = np.load(outcome_path)
        arrival = arrays["arrival_seconds"]
        releases = arrays["release_seconds"]
        if arrival.shape != (48, 1000) or len(cohort) != 48:
            raise ValueError(f"Unexpected particle dimensions for {code}")
        if not np.array_equal(releases, cohort.release_s.to_numpy()):
            raise ValueError(f"Release schedule mismatch for {code}")
        if not np.allclose(np.mean(np.isfinite(arrival), axis=1), cohort.F24):
            raise ValueError(f"F24 mismatch for {code}")
        t50 = np.array([quantile_all_released(row, 0.5) for row in arrival])
        if not np.allclose(t50, cohort.T50_h):
            raise ValueError(f"T50 mismatch for {code}")
        cohort["T90_h"] = [quantile_all_released(row, 0.9) for row in arrival]
        cohort["scenario"] = code
        cohort["E_hinv"] = 1.0 / cohort.T50_h
        peak = reference.loc[reference.Q_ref.idxmax(), "time"]
        cohort["phase"] = np.where(cohort.time <= peak, "Pre-peak", "Post-peak")
        cohort["hours_from_Qpeak"] = (cohort.time - peak).dt.total_seconds() / 3600.0
        matched.append(cohort)

        summaries.append(
            dict(
                scenario=code,
                HCI_peak=float(reference.HCI_ref.max()),
                HCI_mean_hourly=float(reference.HCI_ref.mean()),
                HCI_mean_timeweighted=float(np.trapz(reference.HCI_ref, dx=1) / 168),
                HCI_mean_releases=float(cohort.HCI_ref.mean()),
                RSE=float(reference.RSE_ref.iloc[0]),
                Qpeak_time=str(peak),
                HCI_peak_time=str(reference.loc[reference.HCI_ref.idxmax(), "time"]),
                Cq_ge09_hours=duration_above(reference.Cq_ref.to_numpy(), 0.9),
                HCI_ge05_hours=duration_above(reference.HCI_ref.to_numpy(), 0.5),
                HCI_ge08_hours=duration_above(reference.HCI_ref.to_numpy(), 0.8),
                T50_min=float(cohort.T50_h.min()),
                T50_max=float(cohort.T50_h.max()),
                T50_median=float(cohort.T50_h.median()),
                T50_mean=float(cohort.T50_h.mean()),
                T50_pooled=quantile_all_released(arrival, 0.5),
                T90_pooled=quantile_all_released(arrival, 0.9),
                T90_censored_batches=int(cohort.T90_h.isna().sum()),
                F6_mean=float(cohort.F6.mean()),
                F24_mean=float(cohort.F24.mean()),
                F24_min=float(cohort.F24.min()),
                F24_max=float(cohort.F24.max()),
                rho_HCI_T50=spearman(cohort.HCI_ref, cohort.T50_h),
                rho_HCI_F24=spearman(cohort.HCI_ref, cohort.F24),
                rho_Q_T50=spearman(cohort.Q_ref, cohort.T50_h),
                rho_HCI_Q=spearman(cohort.HCI_ref, cohort.Q_ref),
                rho_diff_HCI_T50=spearman(cohort.HCI_ref.diff().iloc[1:], cohort.T50_h.diff().iloc[1:]),
                T50_min_release=str(cohort.loc[cohort.T50_h.idxmin(), "time"]),
                T50min_minus_Qpeak_h=float((cohort.loc[cohort.T50_h.idxmin(), "time"] - peak).total_seconds() / 3600),
                domain_exit_mean=float(cohort.numerical_domain_exit_fraction.mean()),
                nodata_mean=float(cohort.numerical_nodata_fraction.mean()),
            )
        )

        for tau in np.arange(0, 24.0001, 0.05):
            cdf.append(dict(scenario=code, elapsed_h=tau, F=float(np.mean(arrival <= tau * 3600)), stage="All cohorts"))
        for label, target in (
            ("Pre-peak 24 h", peak - pd.Timedelta(hours=24)),
            ("Near peak", peak),
            ("Post-peak 24 h", peak + pd.Timedelta(hours=24)),
        ):
            index = int(np.argmin(abs(cohort.time - target)))
            representative.append(
                dict(
                    scenario=code,
                    stage=label,
                    batch=index,
                    release_time=str(cohort.time.iloc[index]),
                    offset_h=float((cohort.time.iloc[index] - target).total_seconds() / 3600),
                    T50_h=cohort.T50_h.iloc[index],
                    F24=cohort.F24.iloc[index],
                )
            )
            for tau in np.arange(0, 24.0001, 0.05):
                cdf.append(dict(scenario=code, elapsed_h=tau, F=float(np.mean(arrival[index] <= tau * 3600)), stage=label))
        for source in (cohort_path, outcome_path):
            provenance.append(
                dict(scenario=code, file=source.relative_to(args.particles.parent).as_posix(), sha256=hashlib.sha256(source.read_bytes()).hexdigest())
            )

    summary = pd.DataFrame(summaries)
    all_cohorts = pd.concat(matched, ignore_index=True)
    summary.to_csv(args.output / "summary.csv", index=False)
    all_cohorts.to_csv(args.output / "matched_cohorts.csv", index=False)
    pd.DataFrame(cdf).to_csv(args.output / "cdf.csv", index=False)
    pd.DataFrame(representative).to_csv(args.output / "representative_cohorts.csv", index=False)
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    print(f"Generated validation tables in {args.output}")


if __name__ == "__main__":
    main()
