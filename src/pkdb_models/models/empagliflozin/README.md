# Empagliflozin

Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model of the SGLT2 inhibitor empagliflozin: its oral absorption (a fraction absorbed, food as a slower absorption, the unabsorbed drug excreted with the feces), its distribution into a fast and a deep tissue pool, its glucuronidation in the liver (OATP1B1/1B3 uptake) and in the kidney, the renal excretion of empagliflozin (filtration of the unbound fraction and tubular secretion by OAT3) and of its glucuronides relative to the glomerular filtration rate, and the pharmacodynamics: the inhibition of SGLT2 lowers the renal threshold of glucose (RTG) and increases the urinary glucose excretion (UGE), which depends on the plasma glucose and the GFR. The model describes the effects of the dose, food, hepatic and renal impairment, type 2 diabetes (the fasting plasma glucose) and of the inhibitors of OATP1B1/1B3 and OAT3 (probenecid, rifampicin, gemfibrozil).

The published version of the model is [empagliflozin-model](https://github.com/matthiaskoenig/empagliflozin-model) ([doi:10.5281/zenodo.18480677](https://doi.org/10.5281/zenodo.18480677)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fits, which the commands below run (pkpdsim).
- `models/`: the tissue models (`model_intestine.py`, `model_liver.py`, `model_kidney.py` with the pharmacodynamics) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fits (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the clinical studies (`studies/`), the dose dependency (`misc/`) and the parameter scans of renal function, cirrhosis, dose and food, each again for fasting plasma glucose levels of 5 to 11 mM (`scans/`, with the pharmacokinetic parameters of the studies in `scans/study_data/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`; the base experiment holds the literature based factors of the interactions.
- `fitting/`: the fit problems of the pharmacokinetics (`PK`) and of the urinary glucose excretion on the adopted pharmacokinetics (`PD`), their parameters and the selection of the training, validation (interactions, doses, populations, limitations), outlier and aggregate data.
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics of empagliflozin and its glucuronide and of the UGE.
- `empagliflozin_pk.py`: the pharmacokinetic parameters of the scans, the non-compartmental analysis of `pkpdutils`.
- `tolerances.py`: the absolute tolerances of the integrator per state.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`).
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.empagliflozin --action factory
run_simulation pkdb_models.models.empagliflozin --action list_experiments
run_simulation pkdb_models.models.empagliflozin --action simulate --experiments all
run_simulation pkdb_models.models.empagliflozin --action simulate --experiments all --uncertainty

# fits, adoption of a fit, reports of stored parameters, identifiability
run_fit pkdb_models.models.empagliflozin --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.empagliflozin --adopt PK_<fit id>
run_report pkdb_models.models.empagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.empagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local by default)
python -m pkdb_models.models.empagliflozin.sensitivity.sensitivity_analysis
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fits of the model,
not from numbers in the code. `run_fit pkdb_models.models.empagliflozin --adopt <fit id>`
copies the record of a finished fit (its parameters, problem, profile
likelihood and Fisher information) into `models/output/fits/<fit id>/`, which is
tracked; the ids of the adopted fits are listed in `fits` of the model object
(`__init__.py`). The factory takes the fitted values and their uncertainty from
the records (`models/fits.py`) and writes the uncertainty into the SBML
(`distrib:uncertainty`). `run_simulation ... --uncertainty` runs the
experiments over draws of these distributions and draws the simulations as the
median with a band of quantiles into `results/simulation_uncertainty`. A refit
is: fit, adopt, replace the id in `fits`, rebuild the model. See the
documentation of the repository on fitting (`docs/fitting.md`).
