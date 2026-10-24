"""Pharmacokinetic parameters of the empagliflozin scans.

The parameters are the non-compartmental analysis of `pkpdutils`, named like
the variables of its `NCAResult` (e.g. `auc_last`, `auc_inf_obs`, `cmax`,
`tmax`, `thalf`, `lambda_z`, `cl_f`, `vz_f`), with the unit of every parameter
in `<name>_unit`. Empagliflozin is analysed in plasma with its oral dose, with
the amount excreted in the urine (`Aurine_emp`) and the renal clearance
(`cl_renal`); its glucuronide (EG) and the total of empagliflozin and EG
without a dose, i.e., without clearance and volume. The concentrations below
`LLOQ_FRACTION` of the maximum are below the limit of quantification, so that
the noise of the integrator is not the terminal phase.

A scan has one dimension, the scanned parameter; a run over the draws of the
fitted parameters (`run_simulation pkdb_models.models.empagliflozin
--uncertainty`) adds the dimension of the draws. The parameters are calculated
for every point of every dimension, a row per point with its position along
each dimension (the columns `point` and, with draws, `draw`), see `summarize`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pkpdsim.pk import nca_frame, scan_nca
from pkpdsim.scan import end_values, scan_dimensions, start_values
from pkpdutils import Route, ureg
from sbmlsim import Q
from sbmlsim.result import ScanResult

if TYPE_CHECKING:
    from pkdb_models.models.empagliflozin.experiments.base_experiment import (
        EmpagliflozinSimulationExperiment,
    )

#: oral dose of the scans
DOSE = "PODOSE_emp"
#: the analytes (the names of the column `substance`) -> their concentration
#: and whether the dose is theirs
ANALYTES: dict[str, tuple[str, bool]] = {
    "empagliflozin": ("[Cve_emp]", True),
    "empagliflozin-glucuronide": ("[Cve_eg]", False),
    "total empagliflozin": ("[Cve_emptot]", False),
}
#: cumulative amount of empagliflozin excreted in the urine
URINE = "Aurine_emp"
#: the limit of quantification relative to the maximal concentration
LLOQ_FRACTION = 1e-4


def _point_dimension(res: ScanResult, across: str | None) -> str:
    """Get the dimension of the scan besides the draws `across`.

    Raises:
        ValueError: if the scan has not one dimension besides the draws.
    """
    dims = scan_dimensions(res, across)
    if len(dims) != 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {dims} (draws: {across})."
        )
    return dims[0]


def calculate_empagliflozin_pk(
    experiment: EmpagliflozinSimulationExperiment,
    res: ScanResult,
    across: str | None = None,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of every simulation of a scan.

    The non-compartmental analysis of all simulations at once per analyte
    (`pkpdsim.pk.scan_nca`); empagliflozin with its dose, in the substance
    unit of the concentration which `cl_f` and `vz_f` require, the other
    analytes without a dose. The renal clearance of empagliflozin is the
    amount excreted in the urine until the end of the simulation over
    `auc_last`.

    Args:
        experiment: experiment of the scan, for the molecular weight.
        res: result of the scan, a scan of the oral dose.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        One row per simulation of the scan and analyte, with its position
        (`point`, and `draw` with draws), the points of an analyte in the order
        of the simulations.
    """
    point = _point_dimension(res, across)
    units = res.units
    frames = []
    for substance, (concentration, with_dose) in ANALYTES.items():
        if not with_dose:
            df = nca_frame(
                scan_nca(
                    res,
                    concentration,
                    dose=None,
                    route=None,
                    substance=substance,
                    lloq_fraction=LLOQ_FRACTION,
                )
            )
        else:
            dose_mass = Q(start_values(res, DOSE).values, units[DOSE])
            result = scan_nca(
                res,
                concentration,
                dose=(dose_mass / experiment.Mr.emp).to("mmole"),
                route=Route.ORAL,
                substance=substance,
                lloq_fraction=LLOQ_FRACTION,
            )
            df = nca_frame(result)
            position = tuple(df[dim].to_numpy() for dim in res.dims)
            aurine = ureg.Quantity(
                end_values(res, URINE).values[position], str(units[URINE])
            )
            auc_last = ureg.Quantity(
                df["auc_last"].to_numpy(), df["auc_last_unit"].iloc[0]
            )
            for key, value in [
                (URINE, aurine),
                ("cl_renal", (aurine / auc_last).to("ml/min")),
            ]:
                df[key] = value.magnitude
                df[f"{key}_unit"] = str(value.units)
        df = df.drop(columns="nca_flags")
        df.insert(len(res.dims), "substance", substance)
        frames.append(df)
    names = {point: "point", **({across: "draw"} if across is not None else {})}
    return pd.concat(frames, ignore_index=True).rename(columns=names)


def summarize(
    df: pd.DataFrame, key: str, quantiles: tuple[float, float]
) -> pd.DataFrame:
    """Get the median and the quantiles of a parameter over the draws of every point.

    Args:
        df: the parameters with the columns `point` and `draw`.
        key: the parameter.
        quantiles: the lower and the upper quantile.

    Returns:
        One row per point (the index), the columns `lower`, `median`, `upper`.
    """
    low, high = quantiles
    grouped = df.groupby("point")[key]
    return pd.DataFrame(
        {
            "lower": grouped.quantile(low),
            "median": grouped.median(),
            "upper": grouped.quantile(high),
        }
    ).astype(np.float64)
