Analysis code for the manuscript. All scripts run on institutional compute;
bulk inputs and run outputs are held outside this repository.

Forward model and signal size
  signal_size.py            absorption coefficients and index swing from
                            published optical constants for ice and liquid water
Pilot evidence on field spectra
  index_test.py             tier-1 index on frozen vs unfrozen field spectra
  index_test2.py            group-level permutation tests
  index_tier2.py            weak-band index for saturated organic soils
Freezing-curve reconstruction
  clapeyron_reach2.py       fraction of the curve inside the retention bracket
  sfcc_two_anchor.py        two-anchor reconstruction and error budget
National-scale analyses
  build_retention_dataset.py   md5-gated spectra and retention assembly
  plsr_retention.py            pedon-grouped evaluation, pre-registered
  build_visnir_dataset.py      paired reflectance and property assembly
  applicability.py             predicted detectability across the soil population
  resolution.py                minimum detectable change in liquid fraction
Literature structure
  build_kg.py               modality by soil-state co-occurrence matrix
  analyse_kg.py             temporal trend, property coverage, venue analysis
Multitemporal and field climatology
  get_snotel.py             soil temperature and moisture series near the site
  ladder_design.py          residence-time weighting that sets the laboratory
                            temperature ladder from field climatology
  stac_scan.py              shoulder-season optical scene availability
