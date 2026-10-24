"""The tolerances of the integrator for the empagliflozin model."""

from sbmlsim.model.tolerances import AbsoluteTolerance


def absolute_tolerance(substance: float, other: float) -> AbsoluteTolerance:
    """Get the absolute tolerance of the integrator for the empagliflozin model.

    The species of empagliflozin and its glucuronide are small: the plasma
    concentrations peak at about 1e-5 to 1e-2 mM (0.5 to 800 mg) and fall by
    orders of magnitude in the terminal phase, and the renal threshold for
    glucose depends on the plasma concentration of the kidney with a
    half-maximal effect at 2.5e-6 mM (`KI__RTG_E50`); the amounts in the urine
    and the feces are about 1e-5 to 2 mmole. The glucose in the urine (about
    1e1 to 6e3 mmole) is controlled by the relative tolerance. The other states
    are the doses in mg (`PODOSE_emp`, `IVDOSE_emp`, `cum_dose_emp`, the
    targets of the rate rules).

    Args:
        substance: tolerance of the species, amounts [mmole] and concentrations
            [mM] (multiplied by the volume of the compartment).
        other: tolerance of the targets of the rate rules [mg].

    Returns:
        The absolute tolerance per kind of state.
    """
    return AbsoluteTolerance(amount=substance, concentration=substance, other=other)


#: the absolute tolerance of the simulation experiments
ABSOLUTE_TOLERANCE = absolute_tolerance(substance=1e-14, other=1e-10)
