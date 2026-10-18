# Data sources and access status

## Final HCI calculation

The manuscript figures use the final calculation workbook `三站_20210717零点固定RSE_最终HCI结果.xlsx`. For each station, RSE is fixed at its value at 00:00 on 17 July 2021 and `HCI(t) = RSE_fixed × Cq(t)`.

The normalized tables exported from this workbook are released in `data/processed/hci/`. They contain the hourly discharge used in the calculation, fixed RSE, Cq, HCI, activation parameters, station summaries, and event-stage boundaries. The original Excel workbook is retained outside the repository.

## Hydrological observations

Hourly stage and discharge records were supplied for research use by the Yellow River Conservancy Commission of the Ministry of Water Resources and the Sanmenxia Water Control Project. The local source collections are identified by the authors as `水文资料` and `所需数据——河南大学`. Redistribution permission has not been established, so neither the observations nor hydrology-derived source tables are included here. The exact division of stations and variables between the two providers and the institutional access procedure must be confirmed before publication.

## Sentinel-1A imagery

The following public Sentinel-1A Level-1 IW-mode GRD High Resolution, dual-polarisation (VV/VH), ascending-orbit products were used to delineate the river channel, water extent, and in-channel obstructions:

- `S1A_IW_GRDH_1SDV_20210715T102036_20210715T102101_038787_04939F_4DA0`
- `S1A_IW_GRDH_1SDV_20210715T102101_20210715T102126_038787_04939F_ABCE`
- `S1A_IW_GRDH_1SDV_20210727T102036_20210727T102101_038962_0498E1_3015`
- `S1A_IW_GRDH_1SDV_20210727T102101_20210727T102126_038962_0498E1_9ED3`

The embedded SAFE manifests give acquisition intervals of 10:20:36–10:21:26 UTC on 15 July and 10:20:36–10:21:26 UTC on 27 July 2021, relative orbit 40. The imagery is not duplicated in this repository. It should be obtained from the [Copernicus Data Space Ecosystem Sentinel-1 collection](https://dataspace.copernicus.eu/data-collections/copernicus-sentinel-missions/sentinel-1) under the provider's terms. Earlier draft references to a Google Earth Engine/Landsat maximum-water-extent product do not describe the final data source and have been removed from the availability statement.

## Terrain and roughness

The study-area overview uses a Copernicus DEM product at 30 m spatial resolution. The LISFLOOD-FP terrain combines that DEM with processed high-resolution ground-elevation data supplied by the Henan provincial mapping institute. The provider's formal English name, the high-resolution terrain licence, vertical datum, fusion workflow, and permission to redistribute the fused DEM remain to be confirmed. These terrain files are excluded from this repository.

The model roughness grids were derived by reclassifying a national LUCC dataset. The dataset title, version/year, resolution, access record, licence, and complete class-to-Manning's-n mapping require final confirmation and are not represented as CLCD without source evidence.

## LISFLOOD-FP

The audited simulations used LISFLOOD-FP 5.8.9. The executable is not redistributed. Parameter files are provided in `configs/lisflood/`; the referenced terrain, roughness, initial-condition, and boundary-condition files must be obtained separately with the relevant permissions.
