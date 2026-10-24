# Release notes

## Release notes for empagliflozin-model 0.6.0

Made from pkdb_models 3a0b49144b with pkpdsim 0.11.0.

The changes since 0.5.4.

### Core libraries

pkpdsim 0.11.0: not used by the previous release.

sbmlsim 0.4.0 → 0.13.0 (latest 5 of 428):

- **The experiments of a run run in a pool of worker processes** ([#387](https://github.com/matthiaskoenig/sbmlsim/pull/387)).
- **A figure of the serializer is saved without a second drawing** ([#388](https://github.com/matthiaskoenig/sbmlsim/pull/388)).
- **Formulas are parsed without importing the package petab** ([#389](https://github.com/matthiaskoenig/sbmlsim/pull/389)).
- **The points of a PK observable are analysed in padded groups** ([#390](https://github.com/matthiaskoenig/sbmlsim/pull/390)).
- **The derived units of a model are matched without `areIdentical` of libsbml** ([#383](https://github.com/matthiaskoenig/sbmlsim/pull/383)).

sbmlutils 0.9.6 → 0.17.0 (latest 5 of 260):

- `certifi>=2026.10.19` is a dependency: sbmlutils reaches the network through requests, which verifies with the CA bundle of certifi, and certifi 2026.10.19 removes the distrusted emSign root (GHSA-6q48-6c7w-9gw2) (#510)
- continuous integration, the documentation and the release are built with uv 0.13.4, whose bundled CA store removes the emSign root (GHSA-hq85-5cpj-mr62) (#510)
- model creation returns immutable, grouped construction diagnostics separately from the final validation; the opt-in `strict=True` rejects reported input losses, failed checked libSBML operations and undeclared targets of rules and initial assignments before writing, preserving an existing output file; the default continues to create missing target parameters; `ConstructionError` exposes the findings, annotation resources retained verbatim are informational (#509)
- the SBML parser reports grouped known round-trip losses through `Model.preservation_diagnostics`; factory results carry them independently from the construction and validation findings, also after merging models; the opt-in `strict_preservation=True` of parsing or model creation raises `PreservationError`, and a strict write rejects known source losses before touching its destination; default parsing and writing remain permissive (#509)
- `FactoryResult.validation` exposes the validation result of the final file, after annotations; `create_model(raise_on_error=True)` raises on validation errors while leaving the file available for diagnosis (#509)

### Libraries

| Package | Previous | New |
| --- | --- | --- |
| `libroadrunner` | 2.9.0 | 2.11.0 |
| `matplotlib` | 3.10.8 | 3.11.2 |
| `numpy` | 2.4.3 | 2.5.3 |
| `pandas` | 3.0.1 | 3.0.6 |
| `petab` | 0.8.2 | 0.9.0 |
| `pillow` | 12.1.1 | 12.3.0 |
| `pint` | 0.25.2 | 0.26.1 |
| `pkpdsim` | - | 0.11.0 |
| `pkpdutils` | - | 1.3.0 |
| `pymetadata` | 0.5.11 | 0.8.1 |
| `python-libsbml` | 5.21.1 | 5.21.2 |
| `rich` | 14.3.3 | 15.0.0 |
| `sbmlode` | - | 0.3.0 |
| `sbmlsim` | 0.4.0 | 0.13.0 |
| `sbmlutils` | 0.9.6 | 0.17.0 |
| `scikit-learn` | 1.8.0 | 1.9.1 |
| `scipy` | 1.17.1 | 1.18.1 |
| `typst` | 0.14.8 | 0.15.0 |
| `xarray` | 2026.2.0 | 2026.9.0 |

96 other packages changed.

<details>
<summary>All changed packages</summary>

| Package | Previous | New |
| --- | --- | --- |
| `altair` | 6.0.0 | - |
| `annotated-doc` | 0.0.4 | - |
| `annotated-types` | 0.7.0 | 0.8.0 |
| `antimony` | 3.1.1 | 3.1.3 |
| `anyio` | 4.12.1 | - |
| `asttokens` | 3.0.1 | - |
| `attrs` | 25.4.0 | 26.1.0 |
| `bottleneck` | 1.6.0 | - |
| `certifi` | 2026.2.25 | 2026.10.19 |
| `cffi` | 2.0.0 | - |
| `chardet` | 5.2.0 | 7.6.0 |
| `charset-normalizer` | 3.4.6 | 3.5.2 |
| `click` | 8.3.1 | 8.5.0 |
| `cloudpickle` | - | 3.1.2 |
| `coloredlogs` | 15.0.1 | - |
| `contourpy` | 1.3.3 | 1.4.0 |
| `cryptography` | 46.0.5 | - |
| `decorator` | 5.2.1 | 5.3.1 |
| `deepmerge` | - | 3.0.1 |
| `depinfo` | 2.2.0 | - |
| `elementpath` | 5.1.1 | - |
| `executing` | 2.2.1 | - |
| `fastapi` | 0.135.1 | - |
| `fastobo` | 0.14.1 | - |
| `fonttools` | 4.62.1 | 4.66.1 |
| `google-auth` | 2.49.1 | - |
| `google-auth-oauthlib` | 1.3.0 | - |
| `gspread` | 5.12.4 | - |
| `gspread-pandas` | 3.3.0 | - |
| `h11` | 0.16.0 | - |
| `h5netcdf` | - | 1.8.1 |
| `h5py` | - | 3.16.0 |
| `humanfriendly` | 10.0 | - |
| `idna` | 3.11 | 3.20 |
| `ipython` | 9.11.0 | - |
| `ipython-pygments-lexers` | 1.1.1 | - |
| `jedi` | 0.19.2 | - |
| `joblib` | 1.5.3 | 1.6.0 |
| `kiwisolver` | 1.5.0 | 1.5.1 |
| `latex2mathml` | - | 3.81.1 |
| `libroadrunner` | 2.9.0 | 2.11.0 |
| `lxml` | 6.0.2 | 6.1.3 |
| `markdown` | - | 3.11 |
| `markdown-it-py` | 4.0.0 | 4.2.0 |
| `markupsafe` | 3.0.3 | 3.0.4 |
| `matplotlib` | 3.10.8 | 3.11.2 |
| `matplotlib-inline` | 0.2.1 | - |
| `narwhals` | 2.18.0 | 2.26.0 |
| `networkx` | 3.6.1 | 3.7 |
| `numpy` | 2.4.3 | 2.5.3 |
| `oauthlib` | 3.3.1 | - |
| `packaging` | 26.0 | 26.3 |
| `pandas` | 3.0.1 | 3.0.6 |
| `parso` | 0.8.6 | - |
| `pathspec` | - | 1.1.1 |
| `patsy` | 1.0.2 | - |
| `petab` | 0.8.2 | 0.9.0 |
| `pexpect` | 4.9.0 | - |
| `pillow` | 12.1.1 | 12.3.0 |
| `pint` | 0.25.2 | 0.26.1 |
| `pkdb-analysis` | 0.3.1 | - |
| `pkpdsim` | - | 0.11.0 |
| `pkpdutils` | - | 1.3.0 |
| `platformdirs` | 4.9.4 | 4.12.4 |
| `plotly` | 6.6.0 | - |
| `prompt-toolkit` | 3.0.52 | - |
| `pronto` | 2.7.3 | - |
| `psutil` | 7.2.2 | - |
| `ptyprocess` | 0.7.0 | - |
| `pure-eval` | 0.2.3 | - |
| `pyasn1` | 0.6.3 | - |
| `pyasn1-modules` | 0.4.2 | - |
| `pycparser` | 3.0 | - |
| `pydantic` | 2.12.5 | 2.14.1 |
| `pydantic-core` | 2.41.5 | 2.50.1 |
| `pydoe` | 0.9.7 | - |
| `pygments` | 2.19.2 | 2.21.0 |
| `pymdown-extensions` | - | 12.1 |
| `pymetadata` | 0.5.11 | 0.8.1 |
| `pyparsing` | 3.3.2 | 3.3.3 |
| `pyreadline3` | 3.5.4 | - |
| `python-libnuml` | 1.1.7 | - |
| `python-libsbml` | 5.21.1 | 5.21.2 |
| `python-libsedml` | 2.0.33 | - |
| `python-multipart` | 0.0.22 | - |
| `requests` | 2.32.5 | 2.34.2 |
| `requests-oauthlib` | 2.0.0 | - |
| `rich` | 14.3.3 | 15.0.0 |
| `rpds-py` | 0.30.0 | 2026.9.1 |
| `salib` | 1.5.2 | 1.6.0 |
| `sbmlmath` | - | 0.4.1 |
| `sbmlode` | - | 0.3.0 |
| `sbmlsim` | 0.4.0 | 0.13.0 |
| `sbmlutils` | 0.9.6 | 0.17.0 |
| `scikit-learn` | 1.8.0 | 1.9.1 |
| `scipy` | 1.17.1 | 1.18.1 |
| `setproctitle` | 1.3.7 | - |
| `setuptools` | 82.0.1 | - |
| `stack-data` | 0.6.3 | - |
| `starlette` | 0.52.1 | - |
| `statsmodels` | 0.14.6 | - |
| `threadpoolctl` | 3.6.0 | 3.7.0 |
| `tomli` | - | 2.5.0 |
| `traitlets` | 5.14.3 | - |
| `typing-extensions` | 4.15.0 | 4.16.0 |
| `typing-inspection` | 0.4.2 | 0.4.4 |
| `typst` | 0.14.8 | 0.15.0 |
| `tzdata` | 2025.3 | 2026.5 |
| `urllib3` | 2.6.3 | 2.8.0 |
| `uvicorn` | 0.42.0 | - |
| `wcwidth` | 0.6.0 | - |
| `xarray` | 2026.2.0 | 2026.9.0 |
| `xlsxwriter` | 3.2.9 | - |
| `xmlschema` | 4.3.1 | - |
| `zensical` | - | 0.0.70 |

</details>

### Data

- Removed: `empagliflozin/Macha2013d`
- Changed: `dapagliflozin/Kim2023`, `dapagliflozin/vanderAartvanderBeek2020`, `empagliflozin/Ayoub2017`, `empagliflozin/Brand2012`, `empagliflozin/Chen2015a`, `empagliflozin/Chen2020`, `empagliflozin/Friedrich2013`, `empagliflozin/Hailat2022`, `empagliflozin/Heise2013`, `empagliflozin/Heise2013a`, `empagliflozin/Jiang2023b`, `empagliflozin/Laffel2018`, `empagliflozin/Li2020`, `empagliflozin/Macha2013e`, `empagliflozin/Macha2014`, `empagliflozin/Macha2014b`, `empagliflozin/Macha2014f`, `empagliflozin/Macha2015a`, `empagliflozin/Macha2015b`, `empagliflozin/Sarashina2013`, `empagliflozin/Sarashina2014`, `empagliflozin/Seman2013`, `empagliflozin/Zhao2015`, `hydrochlorothiazide/Heise2015`, `ramipril/Macha2013`

### Analysis

#### Models

- `empagliflozin_body.xml`: added 2 species, 1 reaction, 9 parameters; removed 2 species, 61 parameters; changed values `Kp_emp` 1 → 0.23816, `ftissue_emp` 0.1 → 1.11017
- `empagliflozin_body_flat.xml`: added 2 species, 1 reaction, 19 parameters; removed 13 species, 12 reactions, 77 parameters; changed values `GU__EMPABS_k` 0.01 → 0.0129016, `GU__F_emp_abs` 0.65 → 0.438441, `GU__METEXC_k` 0.001 → 0.000785815, `KI__Mr_glc` 180 → 180.16, `KI__RTG_E50` 7.19e-06 → 4.08786e-06, `KI__RTG_base` 12.5 → 11.7004, `KI__RTG_m_fpg` 0.5 → 2.42163, `KI__RTG_max_inhibition` 0.75 → 0.692919, `Kp_emp` 1 → 0.23816, `ftissue_emp` 0.1 → 1.11017
- `empagliflozin_intestine.xml`: added 2 parameters; removed 7 species, 6 reactions; changed values `EMPABS_k` 0.01 → 0.0129016, `F_emp_abs` 0.65 → 0.438441, `METEXC_k` 0.001 → 0.000785815
- `empagliflozin_kidney.xml`: added 1 species, 9 parameters; removed 2 species, 2 reactions, 8 parameters; changed values `Mr_glc` 180 → 180.16, `RTG_E50` 7.19e-06 → 4.08786e-06, `RTG_base` 12.5 → 11.7004, `RTG_m_fpg` 0.5 → 2.42163, `RTG_max_inhibition` 0.75 → 0.692919
- `empagliflozin_liver.xml`: added 2 parameters; removed 4 species, 4 reactions, 5 parameters

#### Experiments

- Added: `Ayoub2017`, `Brand2012`, `Chen2015a`, `Chen2020`, `DoseDependencyExperiment`, `ElDash2021`, `EmpagliflozinParameterScan`, `Friedrich2013`, `Hailat2022`, `Heise2013`, `Heise2013a`, `Heise2015`, `Jiang2023b`, `Kim2021`, `Kim2023`, `Laffel2018`, `Li2020`, `Macha2013d`, `Macha2013e`, `Macha2014`, `Macha2014b`, `Macha2014f`, `Macha2015a`, `Macha2015b`, `Sarashina2013`, `Sarashina2014`, `Seman2013`, `Zhao2015`, `vanderAartvanderBeek2020`
- Removed: `all`

#### Fits

- `PD`: `PD_20261019_212349__c28ac`, no fit record in the previous release
- `PK`: `PK_20261019_201013__13ddc`, no fit record in the previous release

## Release notes for empagliflozin-model v0.5.4

We are pleased to release a new version of the empagliflozin-model including the 
following changes:

- updates for preprint

Your empagliflozin-model team

## Release notes for empagliflozin-model v0.5.3

We are pleased to release a new version of the empagliflozin-model including the 
following changes:

- model preprint added
- FAIR assessment added
- uv update
- updated github actions

Your empagliflozin-model team

## Release notes for empagliflozin-model v0.5.2

We are pleased to release a new version of the empagliflozin-model including the 
following changes:

- code added

Your empagliflozin-model team

## Release notes for empagliflozin-model v0.5.1

We are pleased to release a new version of the empagliflozin-model including the 
following changes:

- zenodo integration

Your empagliflozin-model team

## Release notes for empagliflozin-model v0.5.0

We are pleased to release a new version of the empagliflozin-model including the 
following changes:

- initial model release

Your empagliflozin-model team

