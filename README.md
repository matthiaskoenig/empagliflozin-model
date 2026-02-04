[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18480677.svg)](https://doi.org/10.5281/zenodo.18480677)
[![GitHub Action](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/python.yml/badge.svg)](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/python.yml)
[![GitHub Action](https://github.com/matthiaskoenig/empagliflozin-model/actions/workflows/docker.yml/badge.svg)](https://github.com/matthiaskoenig/sorafenib-model/actions/workflows/docker.yml)

# empagliflozin model
This repository provides the empagliflozin physiologically based pharmacokinetics/ pharmacodynamics (PBPK/PD) model.

The model is distributed as [SBML](http://sbml.org) format available from [`empagliflozin_body_flat.xml`](./models/empagliflozin_body_flat.xml) with 
corresponding [SBML4humans model report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/empagliflozin-model/main/models/empagliflozin_body_flat.xml) and [model equations](./models/empagliflozin_body_flat.md).

The COMBINE archive is available from [`empagliflozin_model.omex`](./empagliflozin_model.omex).

![model overview](./figures/empagliflozin_model.png)

### Comp submodels
* **liver** submodel [`empagliflozin_liver.xml`](./models/empagliflozin_liver.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/empagliflozin-model/main/models/empagliflozin_liver.xml) and [equations](./models/empagliflozin_liver.md).
* **kidney** submodel [`empagliflozin_kidney.xml`](./models/empagliflozin_kidney.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/empagliflozin-model/main/models/empagliflozin_kidney.xml) and [equations](./models/empagliflozin_kidney.md).
* **intestine** submodel [`empagliflozin_intestine.xml`](./models/empagliflozin_intestine.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/empagliflozin-model/main/models/empagliflozin_intestine.xml) and [equations](./models/empagliflozin_intestine.md).
* **whole-body** submodel [`empagliflozin_body.xml`](./models/empagliflozin_body.xml) with [SBML4humans report](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/empagliflozin-model/main/models/empagliflozin_body.xml) and [equations](./models/empagliflozin_body.md).

## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18480677.svg)](https://doi.org/10.5281/zenodo.18480677)

> Alejandro, J., Elias, M. & König, M. (2026).
> *Physiologically based pharmacokinetic/pharmacodynamic (PBPK) model of empagliflozin.*   
> Zenodo. [https://doi.org/10.5281/zenodo.18480677](https://doi.org/10.5281/zenodo.18480677)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.


## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/empagliflozin-model.git
cd empagliflozin-model
```

#### uv
Run the complete analysis with uv (https://docs.astral.sh/uv/getting-started/installation/):
```bash
uv run run_empagliflozin -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run run_empagliflozin -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/empagliflozin:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_empagliflozin -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```

## Funding
Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B).


© 2024-2026 Jennesse Alejandro, Michelle Elias, and Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
