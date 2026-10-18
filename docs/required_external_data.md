# Required external data

To reproduce the complete workflow, place authorized local copies of the following inputs in a directory outside the Git repository and pass their paths through the command-line interfaces or the reach directories expected by the original particle scripts.

1. Hourly discharge and stage records for HYK, JHT, and GC from 17–24 July 2021.
2. Final hourly HCI source table derived from the paper-final fixed-RSE workbook.
3. LISFLOOD-FP DEM, Manning's-n grid, initial water-depth grid, `.bci` boundary geometry, and `.bdy` hydrograph for each reach.
4. LISFLOOD-FP water-depth and velocity outputs at 600 s intervals.
5. Particle cohort tables and particle-outcome arrays if only the validation summaries and figures are being regenerated.
6. Flood-stage boundary table used to shade the validation figures.

Expected reach directory layouts are documented in the two particle scripts. The public repository intentionally provides no empty files with names that could be mistaken for the restricted originals.

