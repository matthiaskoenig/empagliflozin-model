[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18480677.svg)](https://doi.org/10.5281/zenodo.18480677)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/empagliflozin-model)](https://github.com/matthiaskoenig/empagliflozin-model/releases)
[![Documentation](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/empagliflozin-model/)
[![CI](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/empagliflozin?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/empagliflozin)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of empagliflozin

This repository provides the empagliflozin physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model: the absorption, the tissue distribution, the hepatic and renal glucuronidation and the renal excretion of empagliflozin, and the inhibition of SGLT2, which lowers the renal threshold of glucose and increases the urinary glucose excretion.

Documentation: <https://matthiaskoenig.github.io/empagliflozin-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/empagliflozin-model.git && cd empagliflozin-model
uv sync
uv run run_simulation pkdb_models.models.empagliflozin --action all -r results
```

## How to cite

> Alejandro, J., Elias, M., Babaeva, M. & König, M. (2026). *Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of empagliflozin*. Zenodo. https://doi.org/10.5281/zenodo.18480677

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2024-2026 Jennesse Alejandro, Michelle Elias, Mariia Babaeva and Matthias König
