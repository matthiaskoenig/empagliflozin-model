"""Empagliflozin pharmacokinetics."""

import pandas as pd
from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from sbmlsim.result import XResult
from sbmlutils.log import get_logger


logger = get_logger(__name__)


def calculate_empagliflozin_pk(
    experiment: "EmpagliflozinSimulationExperiment",
    xres: XResult,
) -> pd.DataFrame:
    """
    Calculate empagliflozin parameters.
    Only works for 1D-scans.
    Currently only supporting PO scans.
    """
    Q_ = experiment.Q_

    # scanned dimension
    scandim = xres._redop_dims()[0]

    dose_vec = Q_(xres["PODOSE_emp"].values[0], xres.uinfo["PODOSE_emp"])

    # calculate empagliflozin, empagliflozin-glucuronide, total empagliflozin pharmacokinetic parameters
    pk_dicts = list()
    substance_ids = ["emp", "eg", "emptot"]
    substances = ["empagliflozin", "empagliflozin-glucuronide", "total empagliflozin"]

    t_vec: Q_ = xres.dim_mean("time")
    t_vec = Q_(t_vec.magnitude, xres.uinfo["time"])
    for k_dose, dose in enumerate(dose_vec):

        dose_mmole = dose / experiment.Mr.emp

        for k_sid, sid in enumerate(substance_ids):
            substance = substances[k_sid]
            c_vec = Q_(
                xres[f"[Cve_{sid}]"].sel({scandim: k_dose}).values,
                xres.uinfo[f"[Cve_{sid}]"]
            )
            tcpk = TimecoursePK(
                time=t_vec,
                concentration=c_vec,
                substance=substance,
                dose=dose_mmole,
                ureg=experiment.ureg,
                min_treshold=1e8
            )

            pk_dict = tcpk.pk.to_dict()
            pk_dict["substance"] = substance
            pk_dicts.append(pk_dict)

    return pd.DataFrame(pk_dicts)
