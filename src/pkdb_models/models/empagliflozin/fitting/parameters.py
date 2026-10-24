"""Parameters to optimize in the empagliflozin fits, with their bounds.

The start values are the values of the model, chosen by hand against the
literature (Seman2013, Macha2014f, Chen2015a, Sarashina2013, Heise2013). Not
fitted: the renal fraction of the glucuronidation (`f_ugt_renal`, which the
data of healthy subjects cannot tell from the hepatic one; it was set from the
AUC ratios of hepatic and renal impairment, Macha2014b and Macha2014f, which are
validation data), the dissolution (`GU__Ka_dis_emp`), the fraction unbound
(`KI__fu_emp`), the factors of the interactions (`LI__f_oatp`, `KI__f_oat3`,
`f_ugt`), the hill coefficient (`KI__RTG_gamma`, an exponent below 1 makes the
integration stiff at zero empagliflozin), the splay of the glucose excretion
(`KI__RTG_splay`) and the difference of the mean to the fasting plasma glucose
(`KI__PG_delta`, set by the studies).
"""

from sbmlsim.fit import FitParameter

#: the absorption and the distribution of empagliflozin: tmax 1.0-1.5 hr
#: fasted, food lowers Cmax by 20-29 % and delays tmax by 1.5 hr [Seman2013,
#: Chen2020]; the fraction absorbed about 0.5-0.6 (54.4 % of the radioactivity
#: in the urine, 34.2 % of the dose unchanged in the feces, fecal excretion
#: 37.9 % at 96 hr [Chen2015a]); Vss/F 73.8 l, t1/2 12.4 hr (label), a terminal
#: phase of 16-20 hr in the studies sampled to 96 hr [Macha2014b, Macha2014f]
#: with the deep pool of the rest of the body
parameters_pk_absorption_distribution: list[FitParameter] = [
    FitParameter(
        pid="ftissue_emp",
        start_value=0.45,
        lower_bound=0.01,
        upper_bound=5,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_emp",
        start_value=0.38,
        lower_bound=0.05,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(
        pid="ftissue_deep_emp",
        start_value=0.02,
        lower_bound=1e-3,
        upper_bound=0.5,
        unit="l/min",
    ),
    # the capacity of the deep pool is the volume of the rest of the body
    # (65 l at 75 kg) times Kp_deep: up to 1, i.e. a Vss/F up to about 190 l
    # (73.8 l in the label); a larger pool with a slower exchange (the first
    # fit: Kp_deep 7.4, ftissue_deep 0.004 l/min, profile without an upper
    # bound) is not identified by the data and gives a plateau after 72 hr,
    # 1.5-2 times the concentrations at 96 hr [Macha2014b, Macha2014f]
    FitParameter(
        pid="Kp_deep_emp",
        start_value=0.5,
        lower_bound=0.05,
        upper_bound=1,
        unit="dimensionless",
    ),
    FitParameter(
        pid="GU__EMPABS_k",
        start_value=0.0097,
        lower_bound=1e-3,
        upper_bound=0.1,
        unit="1/min",
    ),
    FitParameter(
        pid="GU__F_emp_abs",
        start_value=0.5,
        lower_bound=0.35,
        upper_bound=0.85,
        unit="dimensionless",
    ),
    FitParameter(
        pid="GU__METEXC_k",
        start_value=0.0012,
        lower_bound=1e-4,
        upper_bound=1e-2,
        unit="1/min",
    ),
    FitParameter(
        pid="GU__FED_f_absorption",
        start_value=0.6,
        lower_bound=0.1,
        upper_bound=1.0,
        unit="dimensionless",
    ),
]
#: the glucuronidation of empagliflozin (hepatic and, `f_ugt_renal` of it,
#: renal) and the renal excretion relative to the GFR: renal clearance of
#: empagliflozin 28-39 ml/min (the filtration of the unbound fraction, 14
#: ml/min, and the secretion by OAT3), the glucuronides about 0.2 of
#: empagliflozin in the plasma and 23 % of the dose in the urine [Chen2015a]
parameters_pk_elimination: list[FitParameter] = [
    FitParameter(
        pid="EMP2EG_k",
        start_value=0.0264,
        lower_bound=1e-3,
        upper_bound=0.5,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__EMPSEC_f",
        start_value=1.2,
        lower_bound=0.05,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__EGEX_f",
        start_value=2.75,
        lower_bound=0.3,
        upper_bound=20,
        unit="dimensionless",
    ),
]
#: the pharmacodynamics, the renal threshold of glucose (RTG): the urinary
#: glucose excretion rises from 1 to 100 mg (healthy 20-79 g/24 hr
#: [Seman2013, Sarashina2013], type 2 diabetes 78-102 g/24 hr with 25-50 mg
#: [Heise2013, Sarashina2014, Macha2014f]); the IC50 of SGLT2 is 3.1 nM
#: (unbound) [Grempler2012], the bounds of the EC50 are 0.1-100 nM (total)
parameters_pd: list[FitParameter] = [
    FitParameter(
        pid="KI__RTG_E50",
        start_value=2.5e-6,
        lower_bound=1e-7,
        upper_bound=1e-4,
        unit="mM",
    ),
    FitParameter(
        pid="KI__RTG_base",
        start_value=10.5,
        lower_bound=8,
        upper_bound=14,
        unit="mM",
    ),
    FitParameter(
        pid="KI__RTG_m_fpg",
        start_value=2.4,
        # positive for the search in log10 space
        lower_bound=0.1,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_max_inhibition",
        start_value=0.69,
        lower_bound=0.3,
        upper_bound=0.95,
        unit="dimensionless",
    ),
]
#: the pharmacokinetics of empagliflozin and EG
parameters_pk: list[FitParameter] = (
    parameters_pk_absorption_distribution + parameters_pk_elimination
)
#: all fitted parameters
parameters_all: list[FitParameter] = parameters_pk + parameters_pd
