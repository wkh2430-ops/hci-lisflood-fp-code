# Required external data

To reproduce the complete workflow, place authorized local copies of the following inputs in a directory outside the Git repository and pass their paths through the command-line interfaces or the reach directories expected by the original particle scripts.

1. Original hourly discharge and stage records for HYK, JHT, and GC from 17–24 July 2021, if the normalized released HCI table must be traced back to provider files.
2. LISFLOOD-FP DEM, Manning's-n grid, initial water-depth grid, `.bci` boundary geometry, and `.bdy` hydrograph for each reach.
3. LISFLOOD-FP water-depth and velocity outputs at 600 s intervals, if particle trajectories must be rerun rather than the released outcomes reanalysed.

The final normalized HCI tables, particle cohort tables, particle-outcome arrays, event-stage boundaries, validation tables, and figure scripts are included. Expected reach layouts for a complete hydrodynamic rerun are documented in the two particle scripts. The public repository intentionally provides no empty files with names that could be mistaken for restricted originals.
