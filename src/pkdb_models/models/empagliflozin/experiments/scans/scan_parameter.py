"""Parameter scans empagliflozin."""

from pathlib import Path
from typing import Dict
import matplotlib.axes
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import pandas as pd
from pint import UnitRegistry
from pint.errors import DimensionalityError
from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot.serialization_matplotlib import FigureMPL, MatplotlibFigureSerializer
from sbmlsim.plot.serialization_matplotlib import plt
from sbmlutils.console import console
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pkdb_models.models.empagliflozin.experiments.base_experiment import EmpagliflozinSimulationExperiment
from pkdb_models.models.empagliflozin.helpers import run_experiments


class EmpagliflozinParameterScan(EmpagliflozinSimulationExperiment):
    """Scan the effect of parameters on pharmacokinetics and pharmacodynamics."""

    study_data_files = {
        "renal_scan": str(Path(__file__).parent / "study_data" / "renal_parameters.tsv"),
        "hepatic_scan": str(Path(__file__).parent / "study_data" / "hepatic_parameters.tsv"),
        "dose_scan": str(Path(__file__).parent / "study_data" / "dose_parameters.tsv"),
    }

    num_points = 15
    glucoses = [5, 6, 7, 8, 9, 10, 11]  # [mM]

    # Substance colors
    substance_colors = {
        "empagliflozin": "black",
        "empagliflozin-glucuronide": "#4A6FA5",
    }

    scan_map = {
        "renal_scan": {
            "parameter": "KI__f_renal_function",
            "default": 1.0,
            "range": np.sort(
                np.append(np.logspace(-1, np.log10(2.0), num=num_points), [1.0])
            ),
            "scale": "log",
            "colormap": "renal_cmap",
            "units": "dimensionless",
            "label": "Renal Function [-]",
        },
        "hepatic_scan": {
            "parameter": "f_cirrhosis",
            "default": 0.0,
            "range": np.linspace(0, 0.9, num=num_points),
            "scale": "linear",
            "colormap": "cirrhosis_cmap",
            "units": "dimensionless",
            "label": "Cirrhosis Degree [-]",
        },
        "dose_scan": {
            "parameter": "PODOSE_emp",
            "range": np.sort(
                np.append(np.logspace(np.log10(0.1), np.log10(800), num=num_points), [100])
            ),
            "scale": "log",
            "colormap": "doses_cmap",
            "units": "mg",
            "label": "Oral Dose [mg]",
        },
    }

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        tcscans = {}

        # Renal simulation
        scan_data = self.scan_map["renal_scan"]
        tcscans["scan_po_renal"] = ScanSim(
            simulation=TimecourseSim(
                Timecourse(
                    start=0,
                    end=24 * 60,  # 24 hours
                    steps=5000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_emp": Q_(50, "mg"),
                    },
                )
            ),
            dimensions=[
                Dimension(
                    "dim_scan",
                    changes={
                        scan_data["parameter"]: Q_(
                            scan_data["range"], scan_data["units"]
                        )
                    },
                ),
            ],
        )

        # Hepatic simulation
        scan_data = self.scan_map["hepatic_scan"]
        tcscans["scan_po_hepatic"] = ScanSim(
            simulation=TimecourseSim(
                Timecourse(
                    start=0,
                    end=24 * 60,  # 24 hours
                    steps=5000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_emp": Q_(50, "mg"),
                    },
                )
            ),
            dimensions=[
                Dimension(
                    "dim_scan",
                    changes={
                        scan_data["parameter"]: Q_(
                            scan_data["range"], scan_data["units"]
                        )
                    },
                ),
            ],
        )

        # Dose simulation
        scan_data = self.scan_map["dose_scan"]
        tcscans["scan_po_dose"] = ScanSim(
            simulation=TimecourseSim(
                Timecourse(
                    start=0,
                    end=24 * 60,  # 24 hours
                    steps=5000,
                    changes={
                        **self.default_changes(),
                    },
                )
            ),
            dimensions=[
                Dimension(
                    "dim_scan",
                    changes={
                        scan_data["parameter"]: Q_(
                            scan_data["range"], scan_data["units"]
                        )
                    },
                ),
            ],
        )

        # scans with different glucose levels for UGE
        for glc in self.glucoses:
            # Renal scan with glucose
            scan_data = self.scan_map["renal_scan"]
            tcscans[f"scan_po_renal_glc{glc}"] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "PODOSE_emp": Q_(50, "mg"),
                            "[KI__fpg]": Q_(glc, "mM"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={
                            scan_data["parameter"]: Q_(
                                scan_data["range"], scan_data["units"]
                            )
                        },
                    ),
                ],
            )

            # Hepatic scan with glucose
            scan_data = self.scan_map["hepatic_scan"]
            tcscans[f"scan_po_hepatic_glc{glc}"] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "PODOSE_emp": Q_(50, "mg"),
                            "[KI__fpg]": Q_(glc, "mM"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={
                            scan_data["parameter"]: Q_(
                                scan_data["range"], scan_data["units"]
                            )
                        },
                    ),
                ],
            )

            # Dose scan with glucose
            scan_data = self.scan_map["dose_scan"]
            tcscans[f"scan_po_dose_glc{glc}"] = ScanSim(
                simulation=TimecourseSim(
                    Timecourse(
                        start=0,
                        end=24 * 60,
                        steps=5000,
                        changes={
                            **self.default_changes(),
                            "[KI__fpg]": Q_(glc, "mM"),
                        },
                    )
                ),
                dimensions=[
                    Dimension(
                        "dim_scan",
                        changes={
                            scan_data["parameter"]: Q_(
                                scan_data["range"], scan_data["units"]
                            )
                        },
                    ),
                ],
            )

        return tcscans

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        """Matplotlib figures."""

        # calculate pharmacokinetic parameters
        self.pk_dfs = self.calculate_empagliflozin_pk()

        return {
            **self.figures_mpl_timecourses(),
            **self.figures_mpl_pharmacokinetics(),
        }

    def _plot_timecourse_scan(self, scan_key: str, sim_key: str, figure_name: str) -> FigureMPL:
        """Create timecourse plot for a parameter scan."""
        sids = [
            "[Cve_emp]",
            "[Cve_eg]",
            "Aurine_emp",
            "Aurine_eg",
            "Afeces_emp",
            "Afeces_eg",
        ]

        scan_data = self.scan_map[scan_key]
        range = scan_data["range"]
        rmin, rmax = range[0], range[-1]

        cmap_name = scan_data["colormap"]
        cmap = getattr(self, cmap_name)

        # 1 row × multiple columns
        nrows = 1
        ncols = len(sids)
        figsize = (6 * ncols, 6 * nrows)

        f, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=figsize,
            sharey=False,
            sharex=True,
            dpi=150,
            squeeze=False,
        )
        plt.subplots_adjust(top=0.88, bottom=0.1, left=0.08, right=0.98, hspace=0.3, wspace=0.35)

        # get data
        Q_ = self.Q_
        xres = self.results[f"task_{sim_key}"]

        # scanned dimension
        scandim = xres._redop_dims()[0]
        parameter_id = scan_data["parameter"]
        par_vec = Q_(
            xres[parameter_id].values[0], xres.uinfo[parameter_id]
        )
        t_vec = xres.dim_mean("time").to(self.units["time"])

        for k_sid, sid in enumerate(sids):
            ax = axes[0, k_sid]

            ymax = 0.0
            for k_par, par in enumerate(par_vec):
                c_vec = Q_(
                    xres[sid].sel({scandim: k_par}).values,
                    xres.uinfo[sid],
                ).to(self.units[sid])

                # update ymax
                cmax = np.nanmax(c_vec.magnitude)
                if cmax > ymax:
                    ymax = cmax

                linewidth = 2.0

                if scan_key == "dose_scan":
                    if scan_data["scale"] == "linear":
                        par_value = par.magnitude if hasattr(par, 'magnitude') else par
                        cvalue = (par_value-rmin)/np.abs(rmax-rmin)
                    elif scan_data["scale"] == "log":
                        par_value = par.magnitude if hasattr(par, 'magnitude') else par
                        cvalue = (np.log10(par_value) - np.log10(rmin)) / np.abs(np.log10(rmax) - np.log10(rmin))
                    color = cmap(cvalue)
                else:
                    par_value = par.magnitude if hasattr(par, 'magnitude') else par
                    default_value = scan_data["default"]

                    if np.isclose(default_value, par_value):
                        color = "black"
                        t_vec_default = t_vec
                        c_vec_default = c_vec
                    else:
                        # calculate color based on position in scan range
                        if scan_data["scale"] == "linear":
                            cvalue = (par_value-rmin)/np.abs(rmax-rmin)
                        elif scan_data["scale"] == "log":
                            cvalue = (np.log10(par_value) - np.log10(rmin)) / np.abs(np.log10(rmax) - np.log10(rmin))
                        color = cmap(cvalue)

                ax.plot(
                    t_vec.magnitude,
                    c_vec.magnitude,
                    color=color,
                    linewidth=linewidth,
                )

            if scan_key != "dose_scan":
                ax.plot(
                    t_vec_default.magnitude,
                    c_vec_default.magnitude,
                    color="black",
                    linewidth=2.0,
                )

            ax: matplotlib.axes.Axes

            ax.set_xlabel(
                f"{self.label_time} [{self.units['time']}]",
                fontdict=self.font,
            )
            ax.tick_params(axis="x", labelsize=self.tick_font_size)
            ax.tick_params(axis="y", labelsize=self.tick_font_size)

            ax.set_ylabel(
                f"{self.labels[sid]} [{self.units[sid]}]",
                fontdict=self.font,
            )

            ax.set_ylim(bottom=0.0, top=1.05 * ymax)
            ax.set_xlim(right=24)

            # Add dose annotation for renal and hepatic scans
            if scan_key in ["renal_scan", "hepatic_scan"]:
                dose_value = Q_(xres["PODOSE_emp"].values[0][0], xres.uinfo["PODOSE_emp"])
                dose_mg = dose_value.to("mg").magnitude

                ax.text(
                    0.98, 0.02,
                    f"Dose: {dose_mg:.0f} mg",
                    transform=ax.transAxes,
                    fontsize=12,
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
                )

        # add colorbar
        cb_ax = f.add_axes(rect=[0.1, 1.0, 0.85, 0.03])
        cb_ax.set_in_layout(True)

        # colorbar range
        if scan_data["scale"] == "linear":
            norm = matplotlib.colors.Normalize(vmin=rmin, vmax=rmax, clip=False)
        elif scan_data["scale"] == "log":
            norm = matplotlib.colors.LogNorm(vmin=rmin, vmax=rmax, clip=False)

        cbar = f.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cb_ax,
            orientation="horizontal",
        )

        # ticks
        ticks = [rmin, rmax]
        if scan_key != "dose_scan" and scan_data["default"] not in ticks:
            ticks.append(scan_data["default"])
            ticks = sorted(ticks)
        cbar.set_ticks(ticks)

        if scan_key == "dose_scan":
            tick_labels = [f"{int(tick)}" for tick in ticks]
        else:
            tick_labels = ticks

        cbar.set_ticklabels(
            tick_labels, **{"size": 15, "weight": "medium"}
        )
        cbar.ax.set_xlabel(
            scan_data["label"], **{"size": 15, "weight": "bold"}
        )

        if scan_key != "dose_scan":
            cbar.ax.axvline(x=scan_data["default"], color="black", linewidth=2)

        return f

    def figures_mpl_timecourses(self) -> Dict[str, FigureMPL]:
        """Timecourse plots for key variables as parameter scan."""
        figures = {}

        # Renal scan
        figures["fig_emp_po__renal__timecourse"] = self._plot_timecourse_scan(
            "renal_scan", "scan_po_renal", "fig_emp_po__renal__timecourse"
        )

        # Hepatic scan
        figures["fig_emp_po__hepatic__timecourse"] = self._plot_timecourse_scan(
            "hepatic_scan", "scan_po_hepatic", "fig_emp_po__hepatic__timecourse"
        )

        # Dose scan
        figures["fig_emp_po__dose__timecourse"] = self._plot_timecourse_scan(
            "dose_scan", "scan_po_dose", "fig_emp_po__dose__timecourse"
        )

        return figures

    def _plot_scan_combined(self, scan_key: str, sim_key: str, pk_parameters: list) -> FigureMPL:
        """Create combined figure with all PK parameters and UGE for a single scan type."""
        scan_data = self.scan_map[scan_key]
        scan_name = scan_key.replace("_scan", "").capitalize()

        ncols = 5
        f, axes = plt.subplots(
            nrows=1, ncols=ncols, figsize=(6 * ncols, 6), dpi=150,
            layout="constrained",
            squeeze=False,
        )
        axes = axes.flatten()

        # Plot PK parameters
        for k, pk_key in enumerate(pk_parameters):
            self._plot_pk_subplot(axes[k], scan_key, sim_key, pk_key, scan_data)

        # Plot UGE
        self._plot_uge_subplot(axes[4], scan_key, sim_key, scan_data)

        return f

    def _molar_mass_for_substance(self, substance: str):
        if substance is None:
            return None
        substance_key = substance.lower()
        return {
            "empagliflozin": self.Mr.emp,
            "empagliflozin-glucuronide": self.Mr.eg,
            "eg": self.Mr.eg,
            "total empagliflozin": self.Mr.emp,
        }.get(substance_key)

    def _convert_pk_quantity(self, y, pk_key: str, substance: str):
        target_unit = self.pk_units[pk_key]
        try:
            return y.to(target_unit)
        except DimensionalityError:
            molar_mass = self._molar_mass_for_substance(substance)
            if molar_mass is None:
                raise
            return (y * molar_mass).to(target_unit)

    def _add_study_data_to_plot(self, ax, scan_key: str, substance: str, parameter: str, study_data_file: str,
                                studies_in_legend: set):
        """Add study data points to PK parameter plots."""
        import pandas as pd

        # Read study data
        df = pd.read_csv(study_data_file, sep='\t', comment='#')

        # Filter for relevant data
        df_filtered = df[
            (df['substance'] == substance) &
            (df['parameter'] == parameter)
            ].copy()

        if len(df_filtered) == 0:
            return None, None

        # Get simulation dose for this scan
        Q_ = self.Q_
        # Derive sim_key from scan_key
        if scan_key == "renal_scan":
            sim_key = "scan_po_renal"
        elif scan_key == "hepatic_scan":
            sim_key = "scan_po_hepatic"
        elif scan_key == "dose_scan":
            sim_key = "scan_po_dose"

        xres = self.results[f"task_{sim_key}"]

        # Get simulation dose (for renal/hepatic scans for scaling)
        if scan_key in ["renal_scan", "hepatic_scan"]:
            sim_dose_value = Q_(xres["PODOSE_emp"].values[0][0], xres.uinfo["PODOSE_emp"]).to("mg").magnitude

        # Get condition to x-axis mapping
        condition_map = None
        if scan_key == "renal_scan":
            # renal_map from base experiment
            condition_map = {
                "normal": self.renal_map["Normal renal function"],
                "mild": self.renal_map["Mild renal impairment"],
                "moderate": self.renal_map["Moderate renal impairment"],
                "severe": self.renal_map["Severe renal impairment"],
            }
        elif scan_key == "hepatic_scan":
            # cirrhosis_map from base experiment
            condition_map = {
                "normal": self.cirrhosis_map["Control"],
                "mild": self.cirrhosis_map["Mild cirrhosis"],
                "moderate": self.cirrhosis_map["Moderate cirrhosis"],
                "severe": self.cirrhosis_map["Severe cirrhosis"],
            }
        elif scan_key == "dose_scan":
            # For dose scan x-position comes directly from the dose value
            pass
        else:
            return None, None

        # Define scalable parameters
        scalable_params = ["cmax", "aucinf"]

        # Get unique studies for different markers
        studies = df_filtered['study'].unique()
        markers = ['s', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'X']
        study_marker_map = {study: markers[i % len(markers)] for i, study in enumerate(studies)}

        # Track maximum y-value for axis scaling
        max_y = 0.0

        # Get substance color
        substance_color = self.substance_colors.get(substance, "black")

        # Plot each data point
        for idx, row in df_filtered.iterrows():
            condition = row['condition'].lower() if pd.notna(row['condition']) else None

            # Get x-position based on scan type
            if scan_key == "dose_scan":
                x_pos = row["dose"]
                dose_ratio = 1.0  # No scaling for dose scan

            else:
                # For renal/hepatic scans, map condition to x-position
                condition_raw = row.get("condition", None)
                condition = condition_raw.lower() if pd.notna(condition_raw) else None

                if condition not in condition_map:
                    continue

                x_pos = condition_map[condition]

                # Get dose and check if scaling is needed
                study_dose = row.get("dose", np.nan)
                dose_ratio = sim_dose_value / study_dose if pd.notna(study_dose) and study_dose > 0 else 1.0

                if pd.notna(study_dose) and study_dose != sim_dose_value:
                    if parameter not in scalable_params:
                        continue

            # Get study marker
            marker = study_marker_map.get(row['study'], 's')

            # Track this study for legend purposes
            if row['study'] not in studies_in_legend:
                studies_in_legend.add(row['study'])

            has_mean = pd.notna(row['mean'])
            has_sd = pd.notna(row['sd'])
            has_median = pd.notna(row['median'])
            has_min = pd.notna(row['min'])
            has_max = pd.notna(row['max'])

            if has_mean:
                # Get value and unit
                y_value = row['mean']
                value_unit = row.get('unit', None)

                # Convert to target units using the same method as simulation data
                if pd.notna(value_unit):
                    y_quantity = Q_(y_value, value_unit)
                    y_quantity = self._convert_pk_quantity(y_quantity, parameter, substance)
                    y_value = y_quantity.magnitude

                # Apply dose scaling if needed
                if parameter in scalable_params and dose_ratio != 1.0:
                    y_value = y_value * dose_ratio

                # Plot mean with optional SD error bars
                if has_sd:
                    y_err = row['sd']
                    if pd.notna(value_unit):
                        y_err_quantity = Q_(y_err, value_unit)
                        y_err_quantity = self._convert_pk_quantity(y_err_quantity, parameter, substance)
                        y_err = y_err_quantity.magnitude

                    if parameter in scalable_params and dose_ratio != 1.0:
                        y_err = y_err * dose_ratio

                    ax.errorbar(x_pos, y_value, yerr=y_err,
                                fmt=marker, markerfacecolor=substance_color, markeredgecolor='black',
                                markeredgewidth=0.8, markersize=8, ecolor='black',
                                capsize=5, capthick=2, alpha=0.9)
                    max_y = max(max_y, y_value + y_err)
                else:
                    ax.plot(x_pos, y_value, marker=marker, markerfacecolor=substance_color,
                            markeredgecolor='black', markeredgewidth=0.8,
                            markersize=8, alpha=0.9, linestyle='')
                    max_y = max(max_y, y_value)

            elif has_median:
                y_value = row['median']
                value_unit = row.get('unit', None)

                # Convert to target units
                if pd.notna(value_unit):
                    y_quantity = Q_(y_value, value_unit)
                    y_quantity = self._convert_pk_quantity(y_quantity, parameter, substance)
                    y_value = y_quantity.magnitude

                # Apply dose scaling if needed
                if parameter in scalable_params and dose_ratio != 1.0:
                    y_value = y_value * dose_ratio

                # Plot median with optional min/max range
                if has_min and has_max:
                    y_min = row['min']
                    y_max = row['max']

                    if pd.notna(value_unit):
                        y_min_quantity = Q_(y_min, value_unit)
                        y_max_quantity = Q_(y_max, value_unit)
                        y_min_quantity = self._convert_pk_quantity(y_min_quantity, parameter, substance)
                        y_max_quantity = self._convert_pk_quantity(y_max_quantity, parameter, substance)
                        y_min = y_min_quantity.magnitude
                        y_max = y_max_quantity.magnitude

                    if parameter in scalable_params and dose_ratio != 1.0:
                        y_min = y_min * dose_ratio
                        y_max = y_max * dose_ratio

                    y_err = [[y_value - y_min], [y_max - y_value]]
                    ax.errorbar(x_pos, y_value, yerr=y_err,
                                fmt=marker, markerfacecolor=substance_color, markeredgecolor='black',
                                markeredgewidth=0.8, markersize=8, ecolor='black',
                                capsize=0, capthick=2, alpha=0.9)
                    max_y = max(max_y, y_max)
                else:
                    ax.plot(x_pos, y_value, marker=marker, markerfacecolor=substance_color,
                            markeredgecolor='black', markeredgewidth=0.8,
                            markersize=8, alpha=0.9, linestyle='')
                    max_y = max(max_y, y_value)

            elif has_min and has_max:
                # Only have range, plot as a bar
                y_min = row['min']
                y_max = row['max']
                value_unit = row.get('unit', None)

                if pd.notna(value_unit):
                    y_min_quantity = Q_(y_min, value_unit)
                    y_max_quantity = Q_(y_max, value_unit)
                    y_min_quantity = self._convert_pk_quantity(y_min_quantity, parameter, substance)
                    y_max_quantity = self._convert_pk_quantity(y_max_quantity, parameter, substance)
                    y_min = y_min_quantity.magnitude
                    y_max = y_max_quantity.magnitude

                if parameter in scalable_params and dose_ratio != 1.0:
                    y_min = y_min * dose_ratio
                    y_max = y_max * dose_ratio

                y_mid = (y_min + y_max) / 2
                y_err = [[y_mid - y_min], [y_max - y_mid]]
                ax.errorbar(x_pos, y_mid, yerr=y_err,
                            fmt=marker, markerfacecolor=substance_color, markeredgecolor='black',
                            markeredgewidth=0.8, markersize=8, ecolor='black',
                            capsize=5, capthick=2, alpha=0.9)
                max_y = max(max_y, y_max)

        return max_y, study_marker_map

    def _add_uge_study_data_to_plot(self, ax, scan_key: str, study_data_file: str, studies_in_legend: set):
        """Add study data points to UGE plots."""
        import pandas as pd

        # Read study data
        df = pd.read_csv(study_data_file, sep='\t', comment='#')

        # Filter for UGE data
        df_filtered = df[df['parameter'] == 'UGE_24'].copy()

        if len(df_filtered) == 0:
            return None, None

        # Get simulation dose for this scan
        Q_ = self.Q_
        # Derive sim_key from scan_key
        if scan_key == "renal_scan":
            sim_key = "scan_po_renal"
        elif scan_key == "hepatic_scan":
            sim_key = "scan_po_hepatic"
        elif scan_key == "dose_scan":
            sim_key = "scan_po_dose"

        # Get simulation dose
        sim_dose_value = None
        if scan_key in ["renal_scan", "hepatic_scan"]:
            glc_first = self.glucoses[0]
            sim_key_glc = f"{sim_key}_glc{glc_first}"
            xres_glc = self.results[f"task_{sim_key_glc}"]
            sim_dose_value = Q_(xres_glc["PODOSE_emp"].values[0][0], xres_glc.uinfo["PODOSE_emp"]).to("mg").magnitude

        # Get condition to x-axis mapping
        condition_map = None
        if scan_key == "renal_scan":
            condition_map = {
                "normal": self.renal_map["Normal renal function"],
                "mild": self.renal_map["Mild renal impairment"],
                "moderate": self.renal_map["Moderate renal impairment"],
                "severe": self.renal_map["Severe renal impairment"],
            }
        elif scan_key == "hepatic_scan":
            condition_map = {
                "normal": self.cirrhosis_map["Control"],
                "mild": self.cirrhosis_map["Mild cirrhosis"],
                "moderate": self.cirrhosis_map["Moderate cirrhosis"],
                "severe": self.cirrhosis_map["Severe cirrhosis"],
            }

        # Map diabetes status to glucose levels and colors
        diabetes_map = {
            "N": 5.0,  # Normal/healthy
            "T2DM": 8.0,  # Type 2 diabetes
            "T1DM": 8.0,  # Type 1 diabetes
        }

        # Get unique studies for different markers
        studies = df_filtered['study'].unique()
        markers = ['s', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'X']
        study_marker_map = {study: markers[i % len(markers)] for i, study in enumerate(studies)}

        max_y = 0.0
        # Plot each data point
        for idx, row in df_filtered.iterrows():
            # Get x-position based on scan type
            if scan_key == "dose_scan":
                x_pos = row['dose']
            else:
                # For renal/hepatic scans, check dose match
                study_dose = row['dose']
                if pd.isna(study_dose) or study_dose != sim_dose_value:
                    continue  # Skip if dose doesn't match

                # Map condition to x-position
                condition_raw = row.get("condition", None)
                condition = condition_raw.lower() if pd.notna(condition_raw) else None

                if condition not in condition_map:
                    continue

                x_pos = condition_map[condition]

            # Get diabetes status and corresponding color
            diabetes_status = row.get('diabetes', 'N')
            if pd.isna(diabetes_status):
                diabetes_status = 'N'

            glucose_level = diabetes_map.get(diabetes_status, 5.0)
            face_color = self.glucose_colors.get(glucose_level, self.glucose_colors[5])

            # Get study marker
            marker = study_marker_map.get(row['study'], 's')

            if row['study'] not in studies_in_legend:
                studies_in_legend.add(row['study'])

            has_mean = pd.notna(row['mean'])
            has_sd = pd.notna(row['sd'])
            has_median = pd.notna(row['median'])
            has_min = pd.notna(row['min'])
            has_max = pd.notna(row['max'])

            if has_mean:
                y_value = row['mean']

                # Plot mean with optional SD error bars
                if has_sd:
                    y_err = row['sd']
                    ax.errorbar(x_pos, y_value, yerr=y_err,
                                fmt=marker, markerfacecolor=face_color, markeredgecolor='black',
                                markeredgewidth=0.8, markersize=8, ecolor='black',
                                capsize=5, capthick=2, alpha=0.9)
                    max_y = max(max_y, y_value + y_err)
                else:
                    ax.plot(x_pos, y_value, marker=marker, markerfacecolor=face_color,
                            markeredgecolor='black', markeredgewidth=0.8,
                            markersize=8, alpha=0.9, linestyle='')
                    max_y = max(max_y, y_value)

            elif has_median:
                y_value = row['median']

                # Plot median with optional min/max range
                if has_min and has_max:
                    y_min = row['min']
                    y_max = row['max']
                    y_err = [[y_value - y_min], [y_max - y_value]]
                    ax.errorbar(x_pos, y_value, yerr=y_err,
                                fmt=marker, markerfacecolor=face_color, markeredgecolor='black',
                                markeredgewidth=0.8, markersize=8, ecolor='black',
                                capsize=0, capthick=2, alpha=0.9)
                    max_y = max(max_y, y_max)
                else:
                    ax.plot(x_pos, y_value, marker=marker, markerfacecolor=face_color,
                            markeredgecolor='black', markeredgewidth=0.8,
                            markersize=8, alpha=0.9, linestyle='')
                    max_y = max(max_y, y_value)

            elif has_min and has_max:
                y_min = row['min']
                y_max = row['max']
                y_mid = (y_min + y_max) / 2
                y_err = [[y_mid - y_min], [y_max - y_mid]]
                ax.errorbar(x_pos, y_mid, yerr=y_err,
                            fmt=marker, markerfacecolor=face_color, markeredgecolor='black',
                            markeredgewidth=0.8, markersize=8, ecolor='black',
                            capsize=5, capthick=2, alpha=0.9)
                max_y = max(max_y, y_max)

        return max_y, study_marker_map

    def _plot_pk_subplot(self, ax, scan_key: str, sim_key: str, pk_key: str, scan_data: dict):
        """Plot a single PK parameter subplot with both substances."""
        Q_ = self.Q_

        # Only show reference line for non-dose scans
        if scan_key != "dose_scan":
            ax.axvline(x=scan_data["default"], color="grey", linestyle="--", linewidth=1.5)

        xres = self.results[f"task_{sim_key}"]
        parameter_id = scan_data["parameter"]
        x_vec = Q_(
            xres[parameter_id].values[0], xres.uinfo[parameter_id]
        )

        substances = ["empagliflozin", "empagliflozin-glucuronide"]
        ymax = 0.0

        # Plot both substances
        for substance in substances:
            df = self.pk_dfs[sim_key]
            df_sub = df[df.substance == substance].copy()

            if len(df_sub) == 0:
                continue

            pk_vec = df_sub[f"{pk_key}"].to_numpy()
            y = Q_(pk_vec, df_sub[f"{pk_key}_unit"].values[0])

            y = self._convert_pk_quantity(y, pk_key, substance)
            color = self.substance_colors[substance]

            substance_label_map = {
                "empagliflozin": "Emp",
                "empagliflozin-glucuronide": "Eg",
            }

            # Plot simulation line
            ax.plot(
                x_vec,
                y,
                marker="o",
                linestyle="-",
                color=color,
                markeredgecolor=color,
                linewidth=2.5,
                markersize=7,
                label=substance_label_map.get(substance, substance),
            )

            ymax = max(ymax, np.nanmax(y.magnitude))

        # Track which studies have been added to legend across all substances
        studies_in_legend = set()
        study_markers = {}  # Map study name to marker symbol

        # Add study data for all substances
        study_max_y = 0.0
        if scan_key in self.study_data_files and self.study_data_files[scan_key] is not None:
            for substance in substances:
                max_y_sub, study_info = self._add_study_data_to_plot(
                    ax, scan_key, substance, pk_key,
                    self.study_data_files[scan_key], studies_in_legend
                )
                if max_y_sub is not None:
                    study_max_y = max(study_max_y, max_y_sub)
                if study_info is not None:
                    study_markers.update(study_info)

        ymax = max(ymax, study_max_y)

        # Formatting
        ax.tick_params(axis="x", labelsize=self.tick_font_size)
        ax.tick_params(axis="y", labelsize=self.tick_font_size)
        ax.set_xlabel(scan_data["label"], fontdict=self.scan_font)
        ax.set_ylabel(
            f"{self.pk_labels[pk_key].capitalize()} [{self.pk_units[pk_key]}]",
            fontdict=self.scan_font,
        )
        ax.set_ylim(bottom=0.0, top=1.05 * ymax)

        if scan_data["scale"] == "log":
            ax.set_xscale("log")
            from matplotlib.ticker import ScalarFormatter
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.get_major_formatter().set_scientific(False)
            if scan_key == "dose_scan":
                ax.set_xticks([0.1, 1, 10, 100, 800])
                ax.set_xticklabels([0.1, 1, 10, 100, 800])

        # Add dose annotation for renal and hepatic scans
        if scan_key in ["renal_scan", "hepatic_scan"]:
            dose_value = Q_(xres["PODOSE_emp"].values[0][0], xres.uinfo["PODOSE_emp"])
            dose_mg = dose_value.to("mg").magnitude
            ax.text(
                0.98, 0.02,
                f"Dose: {dose_mg:.0f} mg",
                transform=ax.transAxes,
                fontsize=12,
                verticalalignment='bottom',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
            )

        # Add legends
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            substance_labels = ["Emp", "Eg"]
            substance_handles = [h for h, l in zip(handles, labels) if l in substance_labels]
            substance_leg_labels = [l for l in labels if l in substance_labels]

            if substance_handles:
                leg1 = ax.legend(
                    substance_handles, substance_leg_labels,
                    title="Substance",
                    loc='upper left',
                    bbox_to_anchor=(0, 0.97),
                    fontsize=10,
                    title_fontsize=11,
                    framealpha=0.9
                )
                ax.add_artist(leg1)

            # Create custom legend handles for studies
            if study_markers:
                import matplotlib.lines as mlines
                study_legend_handles = []
                study_legend_labels = []
                for study_name in sorted(study_markers.keys()):
                    marker = study_markers[study_name]
                    handle = mlines.Line2D(
                        [], [],
                        marker=marker,
                        color='white',
                        markerfacecolor='white',
                        markeredgecolor='black',
                        markeredgewidth=0.8,
                        markersize=8,
                        linestyle='',
                        label=study_name
                    )
                    study_legend_handles.append(handle)
                    study_legend_labels.append(study_name)

                ax.legend(
                    study_legend_handles, study_legend_labels,
                    title="Studies",
                    loc='upper left',
                    bbox_to_anchor=(0.20, 0.97),
                    fontsize=9,
                    title_fontsize=10,
                    framealpha=0.9,
                    ncol=2,
                )

        # Add colorbar
        self._add_colorbar_strip(ax, scan_data)

    def _plot_uge_subplot(self, ax, scan_key: str, sim_key: str, scan_data: dict):
        """Plot UGE subplot with multiple glucose levels."""
        Q_ = self.Q_

        if scan_key != "dose_scan":
            ax.axvline(x=scan_data["default"], color="grey", linestyle="--", linewidth=1.5)

        # Get x-axis values
        sim_key_sample = f"{sim_key}_glc{self.glucoses[0]}"
        xres_sample = self.results[f"task_{sim_key_sample}"]
        parameter_id = scan_data["parameter"]
        x_vec = Q_(
            xres_sample[parameter_id].values[0], xres_sample.uinfo[parameter_id]
        )

        # Plot UGE for each glucose level
        for glc in self.glucoses:
            sim_key_glc = f"{sim_key}_glc{glc}"
            xres_glc = self.results[f"task_{sim_key_glc}"]

            # Get UGE at t=24h
            time_vec = xres_glc.dim_mean("time").magnitude
            t_24h_idx = np.argmin(np.abs(time_vec - 1440))

            # Get UGE values across scan dimension
            scandim = xres_glc._redop_dims()[0]
            uge_values = []
            for k_par in range(len(x_vec)):
                uge_val = Q_(
                    xres_glc["KI__UGE"].sel({scandim: k_par}).values[t_24h_idx],
                    xres_glc.uinfo["KI__UGE"]
                )
                uge_val = uge_val.to(self.units["KI__UGE"])
                uge_values.append(uge_val.magnitude)

            color = self.glucose_colors[glc]
            ax.plot(
                x_vec.magnitude,
                uge_values,
                marker="o",
                linestyle="-",
                color=color,
                linewidth=2,
                markersize=6,
                label=f"{glc} mM"
            )

        # Formatting
        ax.set_xlabel(scan_data["label"], fontdict=self.scan_font)
        ax.set_ylabel(
            f"{self.labels['KI__UGE']} (24h) [{self.units['KI__UGE']}]",
            fontdict=self.scan_font,
        )
        ax.tick_params(axis="x", labelsize=self.tick_font_size)
        ax.tick_params(axis="y", labelsize=self.tick_font_size)
        ymax = ax.get_ylim()[1]

        # Track which studies have been added to legend
        studies_in_legend = set()
        study_markers = {}

        # Add study data if available
        if scan_key in self.study_data_files and self.study_data_files[scan_key] is not None:
            study_max_y, study_info = self._add_uge_study_data_to_plot(
                ax, scan_key, self.study_data_files[scan_key], studies_in_legend
            )
            if study_max_y is not None:
                ymax = max(ymax, study_max_y)
            if study_info is not None:
                study_markers.update(study_info)

        ax.set_ylim(bottom=0.0, top=1.05 * ymax)

        # Set x-scale
        if scan_data["scale"] == "log":
            ax.set_xscale("log")
            from matplotlib.ticker import ScalarFormatter
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.get_major_formatter().set_scientific(False)
            if scan_key == "dose_scan":
                ax.set_xticks([0.1, 1, 10, 100, 800])
                ax.set_xticklabels([0.1, 1, 10, 100, 800])

        # Add dose annotation for renal and hepatic scans
        if scan_key in ["renal_scan", "hepatic_scan"]:
            xres_glc = self.results[f"task_{sim_key}_glc{self.glucoses[0]}"]
            dose_value = Q_(xres_glc["PODOSE_emp"].values[0][0], xres_glc.uinfo["PODOSE_emp"])
            dose_mg = dose_value.to("mg").magnitude
            ax.text(
                0.98, 0.02,
                f"Dose: {dose_mg:.0f} mg",
                transform=ax.transAxes,
                fontsize=12,
                verticalalignment='bottom',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
            )

        # Add legends
        handles, labels = ax.get_legend_handles_labels()
        glucose_handles = [h for h, l in zip(handles, labels) if 'mM' in l]
        glucose_labels = [l for l in labels if 'mM' in l]

        if glucose_handles:
            leg1 = ax.legend(
                glucose_handles, glucose_labels,
                title="Plasma Glucose",
                loc='upper left',
                bbox_to_anchor=(0, 0.97),
                fontsize=9,
                title_fontsize=10,
                framealpha=0.9
            )
            ax.add_artist(leg1)

        if study_markers:
            import matplotlib.lines as mlines
            study_legend_handles = []
            study_legend_labels = []
            for study_name in sorted(study_markers.keys()):
                marker = study_markers[study_name]
                handle = mlines.Line2D(
                    [], [],
                    marker=marker,
                    color='white',
                    markerfacecolor='white',
                    markeredgecolor='black',
                    markeredgewidth=0.8,
                    markersize=8,
                    linestyle='',
                    label=study_name
                )
                study_legend_handles.append(handle)
                study_legend_labels.append(study_name)

            ax.legend(
                study_legend_handles, study_legend_labels,
                title="Studies",
                loc='upper left',
                bbox_to_anchor=(0.24, 0.97),
                fontsize=9,
                title_fontsize=10,
                framealpha=0.9
            )

        # Add colorbar
        self._add_colorbar_strip(ax, scan_data)

    def _add_colorbar_strip(self, ax, scan_data: dict):
        """Add a colorbar strip at the top of a subplot."""
        height_frac = 0.03
        range_vals = scan_data["range"]
        rmin, rmax = range_vals[0], range_vals[-1]
        cmap_name = scan_data["colormap"]
        cmap = getattr(self, cmap_name)

        # Create colorbar normalization
        if scan_data["scale"] == "linear":
            norm = matplotlib.colors.Normalize(vmin=rmin, vmax=rmax, clip=False)
        elif scan_data["scale"] == "log":
            norm = matplotlib.colors.LogNorm(vmin=rmin, vmax=rmax, clip=False)

        # Create inset axes for colorbar strip
        strip_ax = inset_axes(ax, width="100%", height=f"{height_frac * 100:.2f}%", loc="upper center", borderpad=0)

        # Draw colorbar gradient
        if isinstance(norm, matplotlib.colors.LogNorm) or scan_data["scale"] == "log":
            xs = np.geomspace(max(rmin, 1e-12), rmax, 256)
        else:
            xs = np.linspace(rmin, rmax, 256)

        strip_ax.imshow(xs[np.newaxis, :], aspect="auto", cmap=cmap, norm=norm, origin="lower", extent=(0, 1, 0, 1))
        strip_ax.set_axis_off()
        strip_ax.set_zorder(ax.get_zorder() + 1)

    def figures_mpl_pharmacokinetics(self):
        """Visualize dependency of pharmacokinetics parameters."""
        figures = {}

        pk_parameters = ["aucinf", "cmax", "tmax", "thalf"]

        # Scan types
        scan_types = {
            "renal": "renal_scan",
            "hepatic": "hepatic_scan",
            "dose": "dose_scan",
        }

        # Create one figure per scan type
        for scan_label, scan_key in scan_types.items():
            sim_key = f"scan_po_{scan_label}"

            fig_key = f"fig_emp_po__{scan_label}"
            figures[fig_key] = self._plot_scan_combined(
                scan_key, sim_key, pk_parameters
            )

        return figures


if __name__ == "__main__":
    run_experiments(EmpagliflozinParameterScan, output_dir=EmpagliflozinParameterScan.__name__)
