# HCI and LISFLOOD-FP particle-tracking code

This repository contains the analysis code and model-configuration records used to study longitudinal hydrological connectivity in three separate 10 km reaches downstream of the Huayuankou (HYK), Jiahetan (JHT), and Gaocun (GC) gauging stations during the July 2021 lower Yellow River flood.

## Scientific calculation

The final manuscript analysis fixes the relative structural transport efficiency (`RSE`) at its station-specific value at 00:00 on 17 July 2021 and calculates

```text
Cq(t) = 1 / {1 + exp[-beta (Q(t) / Qbank - 1)]}
HCI(t) = RSE_fixed * Cq(t)
```

The fixed RSE values are 0.6469736982 (HYK), 0.9118805156 (JHT), and 0.9244342238 (GC). The exact final hourly outputs used rounded station-specific `k` values in the equivalent expression `Cq(t)=1/{1+exp[-k(Q-Qbank)]}`. Both the manuscript parameterization and the exact reproduction settings are recorded in `configs/stations.json`.

## Repository contents

- `src/hci/calculate_hci.py`: recalculates `Cq(t)` and `HCI(t)` from an authorized hourly discharge table.
- `src/hci/export_final_workbook.py`: exports normalized tables from the paper-final Excel workbook.
- `src/particles/track_hyk_original.py`: HYK passive-particle implementation used in the study.
- `src/particles/track_reach_original.py`: path-configurable JHT/GC passive-particle implementation with wet/dry boundary handling.
- `src/particles/summarize_particles.py`: recalculates cohort `T50`, `T90`, and 6, 12, and 24 h arrival fractions.
- `src/analysis/validation_analysis.py`: pairs release-time HCI with particle outcomes and builds validation source tables.
- `src/analysis/hci_sensitivity.py`: calculates fixed-RSE sensitivity tables for `Qbank`, `beta`, and the activation interval.
- `src/figures/plot_validation.py`: regenerates the validation figure groups from prepared source tables.
- `src/figures/plot_key_figures.py`: compact plotting example for the HCI process and HCI–`T50` relationship.
- `src/lisflood/validate_inputs.py`: audits the packaged LISFLOOD-FP parameter files without executing the model.
- `configs/`: station parameters, particle settings, and the actual LISFLOOD-FP `.par` files.
- `data/processed/hci/`: final hourly HCI, fixed-RSE references, activation parameters, station summaries, and flood-stage boundaries.
- `data/processed/particles/`: 48-cohort particle tables and arrival arrays for the final HYK, JHT, and GC simulations.
- `data/processed/validation/`: matched HCI–particle tables, cumulative-arrival curves, representative cohorts, and three-reach summaries.
- `data/processed/sensitivity/`: fixed-RSE HCI sensitivity tables for `Qbank`, `beta`, and the activation interval.
- `examples/`: synthetic HCI input and expected values.

## Environment

Python 3.10 or later is recommended.

```bash
python -m venv .venv
.venv/Scripts/activate
python -m pip install -r requirements.txt
```

Run the public synthetic example:

```bash
python src/hci/calculate_hci.py \
  --input examples/synthetic_hci_input_and_expected.csv \
  --output results/synthetic_hci_recalculated.csv
```

Rebuild the validation tables, sensitivity tables, and manuscript validation figures from the released processed data:

```bash
python src/analysis/validation_analysis.py
python src/analysis/hci_sensitivity.py
python src/figures/plot_validation.py --output figures
```

The full hydrodynamic simulations require the external files described in `docs/required_external_data.md`. The released processed tables are sufficient to regenerate the HCI and particle-validation figures without rerunning LISFLOOD-FP.

## Reproducibility scope

The repository documents and exposes the HCI calculation, particle-tracking algorithms, processed HCI and particle-result tables, validation analysis, figure code, and LISFLOOD-FP configuration. The normalized final HCI table includes the hourly discharge values used in the published calculation. Original provider workbooks, processed high-resolution terrain, fused model terrain, roughness grids, boundary and initial-condition files, and model flow fields remain excluded. The LISFLOOD-FP executable is also not included.

The four Sentinel-1A IW GRD scenes used to delineate the river channel, water extent, and in-channel obstructions are identified in `docs/data_sources.md`; the imagery should be obtained from the Copernicus Data Space Ecosystem under the original terms.

## Licence

No open-source licence has yet been assigned. See `LICENSE_STATUS.md`. This repository does not grant rights to any third-party dataset or to the LISFLOOD-FP executable.
