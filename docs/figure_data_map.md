# Figure-to-data map

| Figure content | Plotting script | Released input tables |
|---|---|---|
| Hourly HCI process | `src/figures/plot_key_figures.py` | `data/processed/hci/HCI_hourly.csv` |
| HCI and cohort T50 time series | `src/figures/plot_validation.py` | `data/processed/validation/HCI_hourly.csv`, `matched_cohorts.csv`, `summary.csv`, and `data/processed/hci/flood_stage_boundaries.csv` |
| HCI–T50 and HCI–F24 relationships | `src/figures/plot_validation.py` | `data/processed/validation/matched_cohorts.csv`, `summary.csv` |
| Four cumulative-arrival panels | `src/figures/plot_validation.py` | `data/processed/validation/cdf.csv`, `representative_cohorts.csv` |
| Three-reach particle comparison | `src/figures/plot_validation.py` | `data/processed/validation/matched_cohorts.csv`, `summary.csv` |
| Shared-discharge diagnostic | `src/figures/plot_validation.py` | `data/processed/validation/matched_cohorts.csv` |
| Fixed-RSE HCI parameter sensitivity | `src/analysis/hci_sensitivity.py` | `data/processed/sensitivity/HCI_sensitivity_hourly.csv`, `HCI_sensitivity_summary.csv` |

`src/analysis/validation_analysis.py` rebuilds all validation tables from the released hourly HCI table and particle-arrival arrays. The cumulative-arrival denominator is all released particles. The four panels represent all 48 cohorts pooled equally, the nearest release to 24 h before the local discharge peak, the nearest release to the peak, and the nearest release to 24 h after the peak.
