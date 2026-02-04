"""FitParameters for empagliflozin fitting."""
import copy

from sbmlsim.fit import FitParameter

# pharmacokinetic parameters
parameters_pk = [

    # tissue distribution
    FitParameter(
        pid="ftissue_emp",
        lower_bound=0.01,
        start_value=0.1,
        upper_bound=10,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_emp",
        lower_bound=0.1,
        start_value=1,
        upper_bound=10,
        unit="dimensionless",
    ),

    # absorption
    # FitParameter(
    #     pid="GU__F_emp_abs",
    #     lower_bound=0.6,
    #     start_value=0.65,
    #     upper_bound=0.7,
    #     unit="dimensionless",
    # ),
    FitParameter(
        pid="GU__EMPABS_k",
        lower_bound=1E-4,
        start_value=0.01,
        upper_bound=1,
        unit="1/min",
    ),
    # fecal excretion rate
    FitParameter(
        pid="GU__METEXC_k",
        lower_bound=1E-6,
        start_value=0.001,
        upper_bound=0.1,
        unit="1/min",
    ),

    # hepatic metabolism
    # FitParameter(
    #     pid="LI__EMPIM_k",
    #     lower_bound=1,
    #     start_value=10,
    #     upper_bound=100,
    #     unit="1/min",
    # ),
    FitParameter(
        pid="LI__EMP2EG_Vmax",
        lower_bound=1E-3,
        start_value=0.01,
        upper_bound=100,
        unit="mmol/min/l",
    ),
    FitParameter(
        pid="LI__EMP2EG_Km_emp",
        lower_bound=1E-3,
        start_value=0.02,
        upper_bound=1,
        unit="mM",
    ),
    FitParameter(
        pid="LI__EGEX_k",
        lower_bound=1E-3,
        start_value=0.01,
        upper_bound=100,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__EGBIEX_k",
        lower_bound=1E-5,
        start_value=0.01,
        upper_bound=1,
        unit="1/min",
    ),

    # kidney removal
    FitParameter(
        pid="KI__EMPEX_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=10,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__EGEX_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=10,
        unit="1/min",
    ),
]


parameters_pd = [
    FitParameter(
        pid="KI__RTG_E50",
        lower_bound=0.1E-07,
        start_value=2.5e-06,
        upper_bound=50-6,
        unit="mM",
    ),
    FitParameter(
        pid="KI__RTG_base",
        lower_bound=9,
        start_value=12.5,
        upper_bound=14,
        unit="mM",
    ),
    # FitParameter(
    #     pid="KI__RTG_gamma",
    #     lower_bound=1,
    #     start_value=1.01,
    #     upper_bound=4,
    #     unit="dimensionless",
    # ),
    FitParameter(
        pid="KI__RTG_max_inhibition",
        lower_bound=0.2,
        start_value=0.75,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="KI__RTG_m_fpg",
        lower_bound=0.2,
        start_value=1,
        upper_bound=3,
        unit="dimensionless",
    ),

]

parameters_pk_pd = parameters_pk + parameters_pd
