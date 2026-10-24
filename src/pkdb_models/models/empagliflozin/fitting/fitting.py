"""Empagliflozin parameter fitting.

The fit problems of the empagliflozin model: which fit mapping collections
enter a fit and which parameters are adjusted, see `FIT_DEFINITIONS`. The
commands `run_fit`, `run_report` and `run_identifiability` of pkpdsim run them
with the package of the model, a fit is followed by the identifiability
analysis of its parameters (`--no-identifiability` skips it):

    run_fit pkdb_models.models.empagliflozin --subset=PK --runs=3 --cores=3 --seed=1234 --method=DE
    run_fit pkdb_models.models.empagliflozin --subset=PD --runs=3 --cores=3 --seed=1234 --method=DE

The pharmacokinetics are fitted first, the pharmacodynamics on the adopted
values of the pharmacokinetics (`MODEL.fits`). The identifiability analysis is
the profile likelihood of every fitted parameter around the best parameter set
and the Fisher information of the same set, which end up in the
identifiability section of the report of the fit.

The results, figures and the HTML report are written into the fit folder of the
model (`MODEL.fit_dir`), unless `--output_dir` is given. The identifiability of
parameters which were fitted earlier is analysed without fitting again with the
`run_identifiability` command, a report of stored parameters is created with
the `run_report` command:

    run_identifiability pkdb_models.models.empagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json
    run_report pkdb_models.models.empagliflozin --subset=PK results/fit/PK_<fit id>/parameters.json
"""

from pkpdsim.fit import default_fit_settings
from sbmlsim.fit.cli import FitDefinition

from pkdb_models.models.empagliflozin import DATA_PATHS, EMPAGLIFLOZIN_PATH
from pkdb_models.models.empagliflozin.fitting.mapping_collections import (
    f_collections_pd,
    f_collections_pk,
)
from pkdb_models.models.empagliflozin.fitting.parameters import (
    parameters_pd,
    parameters_pk,
)
from pkdb_models.models.empagliflozin.tolerances import absolute_tolerance

#: settings of the fits, stored with the result and read back by the report;
#: the species get a tolerance of their scale, see `tolerances`
FIT_SETTINGS = default_fit_settings(
    absolute_tolerance=absolute_tolerance(substance=1e-12, other=1e-6)
)

#: the fit problems of the empagliflozin model, the keys are the `--subset` argument
FIT_DEFINITIONS: dict[str, FitDefinition] = {
    # the pharmacokinetics of empagliflozin and EG
    "PK": FitDefinition(
        mapping_collections=f_collections_pk,
        parameters=parameters_pk,
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
        settings=FIT_SETTINGS,
    ),
    # the pharmacodynamics: the urinary glucose excretion
    "PD": FitDefinition(
        mapping_collections=f_collections_pd,
        parameters=parameters_pd,
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
        settings=FIT_SETTINGS,
    ),
}
