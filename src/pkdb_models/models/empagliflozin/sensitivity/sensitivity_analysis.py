"""Sensitivity analysis of the empagliflozin model.

The outputs are the pharmacokinetic parameters of empagliflozin and of its
glucuronide after a single oral dose of empagliflozin (the non-compartmental
analysis of `pkpdutils`), the fraction of the dose excreted unchanged in the
urine and the pharmacodynamics: the urinary glucose excretion (UGE) after the
dose and the lowest renal threshold of glucose (RTG). The parameters are the
parameters of the model with bounds of +-15 % around their value. The analyses
are run as a script, the local analysis by default:

    python -m pkdb_models.models.empagliflozin.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.empagliflozin.sensitivity.sensitivity_analysis local sampling

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which one of its groups sets
(`pkpdsim.sensitivity.free_parameters`): `sampling`, over the groups of the
renal function, leaves out `KI__f_renal_function`. The environment variable
`PKPDSIM_SENSITIVITY_N` replaces the numbers of samples of the global analyses
for a quick run, see `pkpdsim.sensitivity.sample_size`.

The model without empagliflozin is stationary (the plasma glucose is constant),
so the simulation of a group starts with the dose. The dose is followed for
`END_HR` hours, long enough for the terminal phase of empagliflozin (a
half-life of 10 to 20 hr) and for the UGE, which returns to its baseline.
"""

from collections.abc import Callable, Mapping

import numpy as np
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from pkpdutils import NCAOptions
from sbmlsim import Q
from sbmlsim.simulation import PK, Custom, Formula, Simulation

from pkdb_models.models.empagliflozin import MODEL, MODEL_PATH
from pkdb_models.models.empagliflozin.experiments.base_experiment import (
    EmpagliflozinSimulationExperiment,
)
from pkdb_models.models.empagliflozin.fitting.parameters import (
    parameters_all as fit_parameters,
)

#: oral dose of empagliflozin of the analysis [mg], the therapeutic dose
DOSE_EMP = 25
#: the hours after the dose
END_HR = 168
#: the limit of quantification of the non-compartmental analysis [mM], about
#: 1e-4 of the maximal plasma concentration after 25 mg (6e-4 mM), so that the
#: noise of the integrator is not the terminal phase
LLOQ = 6e-8
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15

#: outputs of the non-compartmental analysis of empagliflozin -> unit
PK_OUTPUTS_EMP: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "tmax": "min",
    "thalf": "min",
    "vz_f": "l",
    "cl_f": "ml/min",
}
#: outputs of the non-compartmental analysis of the glucuronide, the
#: metabolite of the dose: no clearance or volume -> unit
PK_OUTPUTS_EG: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "thalf": "min",
}
#: the outputs whose unit has the amount of the dose: the dose is a mass
#: (`PODOSE_emp` in mg), so they are divided by the molecular weight (its value
#: at the start, a value per simulation like the parameter of the analysis)
PK_PER_DOSE_MASS: frozenset[str] = frozenset({"vz_f", "cl_f"})


def last(time: np.ndarray, values: dict[str, np.ndarray]) -> float:
    """Get the value of a cumulative amount at the end of the simulation."""
    (x,) = values.values()
    return float(x[time >= 0][-1])


def min_after_dose(time: np.ndarray, values: dict[str, np.ndarray]) -> float:
    """Get the lowest value after the dose."""
    (x,) = values.values()
    return float(np.min(x[time >= 0]))


#: the excretion and the pharmacodynamics after the dose: id -> function,
#: selection, unit
EXTREMES: dict[
    str, tuple[Callable[[np.ndarray, dict[str, np.ndarray]], float], str, str]
] = {
    "emp_aurine": (last, "Aurine_emp", "mmole"),
    "uge": (last, "KI__UGE", "g"),
    "min_rtg": (min_after_dose, "KI__RTG", "mM"),
}

#: the outputs with their ids and names
OUTPUTS: list[SensitivityOutput] = [
    SensitivityOutput("emp_auc_inf_obs", "empagliflozin AUC∞"),
    SensitivityOutput("emp_cmax", "empagliflozin Cmax"),
    SensitivityOutput("emp_tmax", "empagliflozin tmax"),
    SensitivityOutput("emp_thalf", "empagliflozin Half-life"),
    SensitivityOutput("emp_vz_f", "empagliflozin Vz/F"),
    SensitivityOutput("emp_cl_f", "empagliflozin CL/F"),
    SensitivityOutput("emp_aurine", "empagliflozin urine"),
    SensitivityOutput("eg_auc_inf_obs", "glucuronide AUC∞"),
    SensitivityOutput("eg_cmax", "glucuronide Cmax"),
    SensitivityOutput("eg_thalf", "glucuronide Half-life"),
    SensitivityOutput("uge", "UGE"),
    SensitivityOutput("min_rtg", "min RTG"),
]


def _pk_observables(sid: str, outputs: Mapping[str, str]) -> list:
    """Get the non-compartmental analysis of a substance in plasma and its outputs."""
    pk_id = f"pk_{sid}"
    return [
        PK(
            pk_id,
            f"[Cve_{sid}]",
            dose="PODOSE_emp",
            route="oral",
            options=NCAOptions(lloq=LLOQ),
            parameters=list(outputs),
        ),
        *(
            Formula(
                f"{sid}_{key}",
                f"{pk_id}.{key} / at(Mr_emp, 0)"
                if key in PK_PER_DOSE_MASS
                else f"{pk_id}.{key}",
                unit=unit,
            )
            for key, unit in outputs.items()
        ),
    ]


#: the observables of a simulation, the outputs with their units; the ids of
#: the analyses (`pk_emp`, `pk_eg`) are no species of the model
OBSERVABLES = [
    *_pk_observables("emp", PK_OUTPUTS_EMP),
    *_pk_observables("eg", PK_OUTPUTS_EG),
    *(
        Custom(key, function, unit, symbols=[selection])
        for key, (function, selection, unit) in EXTREMES.items()
    ),
]

#: the output times [hr]: every 0.02 hr in the first day after the dose
#: (absorption and maximum), every 0.5 hr after it
TIMES_HR: tuple[float, ...] = tuple(
    np.unique(
        np.concatenate(
            [
                np.linspace(0, 24, 1201),
                np.linspace(24, END_HR, 2 * (END_HR - 24) + 1),
            ]
        )
    ).tolist()
)


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: a single oral dose of empagliflozin at time 0."""
    return Simulation(
        start=0,
        end=END_HR,
        times=TIMES_HR,
        time_unit="hr",
        preinit_changes={**changes, "PODOSE_emp": Q(DOSE_EMP, "mg")},
    )


#: the groups the analyses are run for: the control and the stages of renal
#: impairment of the simulation experiments
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {
                    "KI__f_renal_function": (
                        EmpagliflozinSimulationExperiment.renal_map[name]
                    )
                }
            ),
            EmpagliflozinSimulationExperiment.renal_colors[name],
        )
        for uid, name in [
            ("mildRI", "Mild renal impairment"),
            ("modRI", "Moderate renal impairment"),
            ("sevRI", "Severe renal impairment"),
        ]
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds and references."""
    parameters = parameters_from_sbml(
        MODEL_PATH,
        exclude_ids={
            # molecular weights and conversion factors
            "Mr_emp",
            "Mr_eg",
            "KI__Mr_glc",
            "KI__cf_mg_per_g",
            # reference values
            "BW_ref",
            # unchangeable values
            "FQlu",
            "FVhv",
            "FVpo",
            # dosing parameters
            "PODOSE_emp",
            "IVDOSE_emp",
            "ti_emp",
            "Ri_emp",
            # unused volumes
            "Vurine",
            "Vfeces",
            "Vstomach",
            # the switch of food and its effect, without effect in the fasted
            # state of the analysis
            "GU__fed",
            "GU__FED_f_absorption",
        },
        exclude_na=True,
        exclude_zero=True,
    )

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters,
        [(fp.pid, np.nan, np.nan, ParameterType.FIT) for fp in fit_parameters],
    )
    # the scaling factor of the renal function
    set_bounds(
        parameters,
        [
            (
                "KI__f_renal_function",
                1 - BOUNDS_FRACTION,
                1 + BOUNDS_FRACTION,
                ParameterType.SCALING,
            )
        ],
    )

    # the references of the values from the literature
    references = {
        "BW": r"\cite{Ogden2004, Jones2013, Thompson2009, Brown1997}",
        "FQgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FQh": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FQki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVar": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVli": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FVlu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVve": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
        "KI__f_renal_function": r"\cite{Stevens2024}",
    }
    by_uid = {p.uid: p for p in parameters}
    for uid, reference in references.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script, with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS[:1], delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "fast": FastAnalysis(groups=SENSITIVITY_GROUPS[:1], n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local",),
)


if __name__ == "__main__":
    SENSITIVITY.main()
