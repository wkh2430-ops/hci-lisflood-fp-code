# Processed-data dictionary

## HCI tables

- `HCI_hourly.csv`: `time`, station code, discharge (`Q_m3_s`), fixed RSE, activation factor Cq, and HCI. Times are local source times; the source files contain no UTC offset.
- `fixed_RSE_reference.csv`: station-specific reference values at 00:00 on 17 July 2021.
- `Cq_parameters.csv`: bankfull discharge, lower and upper activation-reference discharges, and the rounded legacy `k` used in the final workbook.
- `HCI_station_summary.csv`: event minima, maxima, means, and peak timing.
- `flood_stage_boundaries.csv`: discharge-based rising, peak, and recession boundaries used for figure shading.

## Particle tables

- `particle_cohorts.csv`: one row per release cohort. `T50_h` uses all 1000 released particles as the denominator. `F6`, `F12`, and `F24` are the fractions of all released particles arriving within 6, 12, and 24 h.
- `particle_outcomes.npz`: 48 × 1000 arrival times in seconds, particle states, release times, and retained final coordinates needed to recompute cumulative-arrival curves and quantiles.
- `matched_cohorts.csv`: release-time HCI joined to cohort outcomes for HYK, JHT, and the final GC simulation.
- `cdf.csv`: cumulative arrival fraction `F` from 0 to 24 h at 0.05 h intervals. `stage` identifies the pooled and three representative-cohort panels.
- `representative_cohorts.csv`: exact cohort index, release time, and time offset used for each representative panel.
- `summary.csv`: event and cohort statistics, including pooled T50/T90, cohort T50, F24, and descriptive Spearman coefficients.

## Sensitivity tables

- `HCI_sensitivity_hourly.csv`: fixed-RSE hourly scenarios for Qbank factors, beta values, and activation intervals.
- `HCI_sensitivity_summary.csv`: scenario means and maxima, with mean-HCI change relative to the exact final-workbook baseline.

Missing values indicate unavailable or censored quantities; they are not converted to zero.
