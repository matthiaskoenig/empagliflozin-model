"""Parameter fit problems for empagliflozin.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the training data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The training data are the oral doses of healthy subjects and of patients with
type 2 diabetes and a normal renal function, fasted, fed (the effect of food on
the absorption is fitted) or fasting not reported, without a coadministered
drug, single doses and exactly encoded repeated doses. The validation data has
four parts, each a different question: the interactions
(`INTERACTION_VALIDATION`, the arms with a coadministered drug, simulated with
the literature based factors of OAT3 and OATP1B1 or without a factor where the
studies found no interaction), the doses (`DOSE_VALIDATION`, the doses above
the range of the training), the populations (`POPULATION_VALIDATION`, the
impairment of the liver and the kidneys, East Asian patients and adolescents
with type 2 diabetes) and the limitations (`LIMITATION_VALIDATION`, data of a
dose, a meal or a population the model does not describe). The summaries of
data which is fitted one by one, and the same subjects counted twice, are
excluded (`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import mapping_metadata, observable_filter, observable_id
from pkpdsim.metadata import Coadministration, Fasting, Health, Route
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.empagliflozin import DATA_PATHS, EMPAGLIFLOZIN_PATH
from pkdb_models.models.empagliflozin.experiments.studies import (
    Ayoub2017,
    Brand2012,
    Chen2015a,
    Chen2020,
    ElDash2021,
    Friedrich2013,
    Hailat2022,
    Heise2013,
    Heise2013a,
    Heise2015,
    Jiang2023b,
    Kim2021,
    Kim2023,
    Laffel2018,
    Li2020,
    Macha2013d,
    Macha2013e,
    Macha2014,
    Macha2014b,
    Macha2014f,
    Macha2015a,
    Macha2015b,
    Sarashina2013,
    Sarashina2014,
    Seman2013,
    Zhao2015,
    vanderAartvanderBeek2020,
)

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share (Brand2012 and Friedrich2013,
#: Heise2013 and Heise2013a), so a named mapping is named with its study
type NamedMappings = dict[str, set[str]]

EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
    Ayoub2017,
    Brand2012,
    Chen2015a,
    Chen2020,
    ElDash2021,
    Friedrich2013,
    Hailat2022,
    Heise2013,
    Heise2013a,
    Heise2015,
    Jiang2023b,
    Kim2021,
    Kim2023,
    Laffel2018,
    Li2020,
    Macha2013d,  # the empagliflozin groups of the ramipril study Macha2013
    Macha2013e,
    Macha2014,
    Macha2014b,
    Macha2014f,
    Macha2015a,
    Macha2015b,
    Sarashina2013,
    Sarashina2014,
    Seman2013,
    vanderAartvanderBeek2020,
    Zhao2015,
]

#: observables of the pharmacokinetics
PK_OBSERVABLES: set[str] = {
    "Cve_emp",
    "Cve_eg",
    "Cve_emptot",
    "Aurine_emp",
    "Aurine_eg",
    "Aurine_emptot",
    "Afeces_emp",
    "Afeces_emptot",
}

#: health of the training data: healthy subjects and patients with type 2
#: diabetes (or hypertension) with a normal renal function
TRAINING_HEALTH: frozenset[Health] = frozenset(
    {Health.HEALTHY, Health.T2DM, Health.HYPERTENSION}
)

#: fasting of the training data, fasted, fed or not reported: the effect of
#: food on the absorption (`GU__FED_f_absorption`) is fitted with the fed arms
TRAINING_FASTING: frozenset[Fasting] = frozenset(
    {Fasting.FASTED, Fasting.FED, Fasting.NR}
)

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {
    # the glucuronides as the difference of the plasma radioactivity and the
    # empagliflozin of LC-MS/MS: zero up to 1 hr and 2-2.2 times below the
    # glucuronides of Table 2, which are fitted (`fm_po_emp_EG_plasma_Tab2`)
    "Chen2015a": {"fm_po_emp_EG_plasma"},
    # the urinary glucose excretion: the curated values are 5-9 times below
    # the ones of Table 2, which is not consistent in itself (g/dl, cumulative
    # values which fall)
    "ElDash2021": {
        "po_EMP10_uge_morning",
        "po_EMP10_uge_evening",
        *(
            f"po_EMP10_uge_{intervention}_{subject}"
            for intervention in ElDash2021.interventions
            for subject in ElDash2021.subjects
        ),
    },
    # tmax 4.8 hr with a second peak, t1/2 4.0-4.5 hr, AUC 41 % of Chen2020,
    # the fed exposure 3 times below the fasted one against the text, the
    # units of Fig3 unclear, the fasting dose with 48 g of glucose
    "Hailat2022": {
        "fm_hailat2022_fasting_test",
        "fm_hailat2022_fasting_reference",
        "fm_hailat2022_fed_test",
        "fm_hailat2022_fed_reference",
    },
    # the 24 hr values after a single dose of 10 mg are above the 4-6 hr values
    "vanderAartvanderBeek2020": {"fm_emp10_S1", "fm_emp10_S6"},
}

#: the validation of the interactions: the arms with a coadministered drug.
#: Probenecid (OAT3), rifampicin (OATP1B1) and gemfibrozil (OATP1B1, OAT3) are
#: simulated with literature based factors which are not fitted
#: (`interaction_changes` of the experiments), the other drugs without a
#: factor, the studies found no interaction. `filter_training` lets the
#: validation data through for this reason
INTERACTION_VALIDATION: NamedMappings = {
    # sitagliptin, AUCtau,ss ratio 1.10
    "Brand2012": {"fm_empagliflozin_EMP50_D10", "fm_empagliflozin_urine_EMP50_D10"},
    # linagliptin, AUCtau,ss ratio 1.02
    "Friedrich2013": {"fm_empagliflozin_EMP50_D7", "fm_empagliflozin_urine_EMP50_D7"},
    # hydrochlorothiazide and torasemide
    "Heise2015": {"fm_po_emp25_emp25_hctz25", "fm_po_emp25_emp25_tor5"},
    # lobeglitazone
    "Kim2021": {"fm_po_emp25_E__L"},
    # evogliptin
    "Kim2023": {"fm_emp25_EV5__EP25", "fm_uge_EV5__EP25"},
    # verapamil and ramipril, AUC ratios 1.03 and 0.97
    "Macha2013d": {
        "fm_po_emp25_emp25_single__ver120_single",
        "fm_po_emp25_emp25_multi__ram5_multi",
    },
    # warfarin, AUC ratio 1.01
    "Macha2013e": {"fm_po_emp25_emp__war"},
    # gemfibrozil, rifampicin and probenecid with their factors, AUC ratios
    # 1.59, 1.35 and 1.53
    "Macha2014": {
        "fm_po_EMP25_GEM",
        "fm_po_EMP10_RIF",
        "fm_po_EMP10_PROB",
        "fm_po_EMP10_PROB_urine",
    },
    # pioglitazone, AUC ratio 1.00
    "Macha2015b": {"fm_po_emp50_emp50__piog45"},
}

#: the validation of the doses: the doses above the range of the training,
#: whose urinary glucose excretion is the one of a saturated inhibition of
#: SGLT2
DOSE_VALIDATION: NamedMappings = {
    "Seman2013": {"fm_uge_EMP400", "fm_uge_EMP800"},
}

#: the validation of the populations: the impairment of the liver (the shared
#: cirrhosis model) and the kidneys (the renal function of the GFR), East
#: Asian and adolescent patients with type 2 diabetes, whose exposure the model
#: describes with the body weight and the renal function alone
POPULATION_VALIDATION: NamedMappings = {
    # adolescents with type 2 diabetes, BW 90-111 kg, eGFR 157-178 ml/min
    "Laffel2018": {
        "fm_laffel2018_po_EMP5",
        "fm_laffel2018_po_EMP10",
        "fm_laffel2018_po_EMP25",
    },
    # hepatic impairment, Child-Pugh A, B and C
    "Macha2014b": {
        f"fm_{observable}_{group}"
        for observable in [
            "empagliflozin",
            "empagliflozin_urine",
            "uge",
            "eg_urine",
            "emptot_urine",
        ]
        for group in ["mild", "moderate", "severe"]
    },
    # renal impairment of patients with type 2 diabetes
    "Macha2014f": {
        f"fm_{observable}_{group}"
        for observable in ["empagliflozin", "empagliflozin_urine", "uge"]
        for group in ["mild", "moderate", "severe"]
    },
    # Japanese patients with type 2 diabetes, normal and impaired renal function
    "Sarashina2014": {
        f"fm_{observable}_{group}"
        for observable in ["empagliflozin", "empagliflozin_urine", "uge"]
        for group in ["normal", "mild", "moderate", "severe"]
    },
    # patients with type 2 diabetes and albuminuria, individual data
    "vanderAartvanderBeek2020": {
        "fm_emp10_S2",
        "fm_emp10_S3",
        "fm_emp10_S4",
        "fm_emp10_S5",
    },
    # Chinese patients with type 2 diabetes, day 1 and day 9
    "Zhao2015": {
        f"fm_empagliflozin_{observable}{dose}"
        for observable in ["SING", "MULTI", "urine_SING", "urine_MULTI"]
        for dose in ["10", "25"]
    },
}

#: the validation of the limitations of the model: valid data of a dose, a meal
#: or a population the model does not describe, which would distort the fitted
#: parameters of the processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    # the oral solution: 1.6 times the exposure of the tablets (AUC 15600 vs
    # about 9900 nM*h after 50 mg) and fe 0.29 (tablets 0.16-0.20), which one
    # fraction absorbed for every formulation cannot give (the first fit
    # predicted 0.64 of the plasma and 0.72 of the urine); the feces and the
    # glucuronides of the mass balance stay training data
    "Chen2015a": {
        "fm_po_emp_EMP_plasma",
        "fm_po_emp_EMPTOT_plasma",
        "fm_po_emp_EMP_urine",
        "fm_po_emp_EMPTOT_urine",
    },
    # the evening dose with dinner: the model has no time of day, the same
    # subjects as the morning arm
    "ElDash2021": {"po_EMP10_empagliflozin_evening"},
    # the fasting arm with about 100 g of oral glucose with and after the
    # dose, tmax 4 hr like the fed arm
    "Li2020": {"fm_po_emp25_ref_fasting"},
    # end stage renal disease, the dialysis is not modelled
    "Macha2014f": {
        "fm_empagliflozin_esrd",
        "fm_empagliflozin_urine_esrd",
        "fm_uge_esrd",
    },
    # 400 and 800 mg, the exposure is less than dose proportional, the
    # absorption of the model is linear
    "Seman2013": {
        "fm_empagliflozin_EMP400",
        "fm_empagliflozin_EMP800",
        "fm_empagliflozin_urine_EMP400",
        "fm_empagliflozin_urine_EMP800",
    },
}

#: the summaries of data which the fits use one by one, and the same subjects
#: counted twice, excluded from every fit
AGGREGATE_MAPPINGS: NamedMappings = {
    # the test formulation, same subjects as the reference, bioequivalent
    "Chen2020": {"fm_po_emp10_Fasting_T", "fm_po_emp10_Fed_T"},
    # the individual data of the subjects, the means are fitted
    "ElDash2021": {
        f"po_EMP10_empagliflozin_{intervention}_{subject}"
        for intervention in ElDash2021.interventions
        for subject in ElDash2021.subjects
    },
    # the cocrystal test formulation, same subjects as the reference
    "Jiang2023b": {"fm_emp25_test"},
    # the 10 mg arm of the probenecid study drawn a second time (Figure 2C)
    "Macha2014": {"fm_po_EMP10_P"},
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(
            mappings.get(study, set())
            for mappings in [
                INTERACTION_VALIDATION,
                DOSE_VALIDATION,
                POPULATION_VALIDATION,
                LIMITATION_VALIDATION,
            ]
        )
    )
    for study in INTERACTION_VALIDATION
    | DOSE_VALIDATION
    | POPULATION_VALIDATION
    | LIMITATION_VALIDATION
}


def study(fit_mapping: FitMapping) -> str:
    """Get the study of a fit mapping, the id of its experiment."""
    return fit_mapping.reference.experiment.sid


def named(
    mappings: NamedMappings, fit_mapping_key: str, fit_mapping: FitMapping
) -> bool:
    """Check whether a fit mapping is one of the named mappings of its study."""
    return fit_mapping_key in mappings.get(study(fit_mapping), set())


def check_named(mappings: NamedMappings, fit_mappings: FitMappings) -> None:
    """Check that the named mappings are fit mappings of their studies.

    Raises:
        ValueError: for a study which is no experiment, or a key which is no fit
            mapping of its study.
    """
    for sid, keys in mappings.items():
        experiment = fit_mappings.experiments.get(sid)
        if experiment is None:
            raise ValueError(f"'{sid}' is no experiment of the fits.")
        if unknown := sorted(keys - set(experiment._fit_mappings)):
            raise ValueError(f"'{unknown}' are no fit mappings of '{sid}'.")


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No summaries of data which is fitted one by one, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


def filter_training(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only the training data, besides the validation data.

    Oral empagliflozin without a coadministered drug in healthy subjects and
    patients with type 2 diabetes and a normal renal function, fasted, fed or
    fasting not reported.
    """
    if named(VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping):
        return True
    metadata = mapping_metadata(fit_mapping)
    return (
        metadata.route == Route.PO
        and metadata.coadministration == Coadministration.NONE
        and metadata.health in TRAINING_HEALTH
        and metadata.fasting in TRAINING_FASTING
    )


#: Only the pharmacokinetic data.
filter_pk = observable_filter(PK_OBSERVABLES)


def filter_pd(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only the pharmacodynamic data, the urinary glucose excretion."""
    return observable_id(fit_mapping) not in PK_OBSERVABLES


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the training data, the mappings which fail them are
    excluded, and so are the summaries of data which is fitted one by one
    (`AGGREGATE_MAPPINGS`) and the data which is not training data
    (`filter_training`). The outliers are `OUTLIER_MAPPINGS`, the validation
    data is `VALIDATION_MAPPINGS`. Neither kind enters the cost; both are
    simulated and evaluated with the training data, so a report says how the
    fitted parameters describe data they were not fitted on.

    Args:
        filters: filters of the training data.
        validation: the keys of the validation data by study.

    Returns:
        The fit mapping collections of all kinds by experiment id.

    Raises:
        ValueError: if an outlier is a mapping of no study or of several.
    """
    fit_mappings = FitMappings(
        experiment_classes=EXPERIMENT_CLASSES,
        base_path=EMPAGLIFLOZIN_PATH,
        data_path=DATA_PATHS,
    )
    for mappings in [OUTLIER_MAPPINGS, AGGREGATE_MAPPINGS, validation]:
        check_named(mappings, fit_mappings)
    # sbmlsim selects the outliers by their keys alone, which must name one
    # mapping of all studies
    outliers = {key for keys in OUTLIER_MAPPINGS.values() for key in keys}
    for key in outliers:
        studies = [
            sid
            for sid, experiment in fit_mappings.experiments.items()
            if key in experiment._fit_mappings
        ]
        if len(studies) != 1:
            raise ValueError(f"The outlier '{key}' is a mapping of {studies}.")
    selection = fit_mappings.select(
        filters=[
            *(filters if isinstance(filters, list) else [filters]),
            filter_aggregates,
            filter_training,
        ],
        outliers=outliers,
        validation=lambda key, fit_mapping: named(validation, key, fit_mapping),
        print_info=False,
    )
    # the counts per experiment and kind, not the single mappings
    selection.print(detail=False)
    # the validation mappings by the question they answer
    selected = {
        (sid, key)
        for sid, kinds in selection.kinds.items()
        for key, kind in kinds.items()
        if kind is MappingKind.VALIDATION
    }

    def pairs(mappings: NamedMappings) -> set[tuple[str, str]]:
        return {(sid, key) for sid, keys in mappings.items() for key in keys}

    console.print(
        "validation: "
        + ", ".join(
            f"{len(selected & pairs(mappings))} {question}"
            for question, mappings in [
                ("interactions", INTERACTION_VALIDATION),
                ("doses", DOSE_VALIDATION),
                ("populations", POPULATION_VALIDATION),
                ("limitations", LIMITATION_VALIDATION),
            ]
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """Pharmacokinetic data."""
    return mapping_collections([filter_pk])


def f_collections_pd() -> dict[str, list[FitMappingCollection]]:
    """Pharmacodynamic data, the urinary glucose excretion."""
    return mapping_collections([filter_pd])


if __name__ == "__main__":
    for f in [f_collections_pk, f_collections_pd]:
        console.rule(title=f.__name__, align="left", style="white")
        f()
