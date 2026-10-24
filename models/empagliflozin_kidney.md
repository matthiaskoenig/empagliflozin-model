# Model for the renal excretion of empagliflozin and EG and the urinary glucose excretion.

Model `empagliflozin_kidney`, SBML Level 3 Version 2, read from empagliflozin\_kidney.xml, written by sbmlode 0.3.0.

Model for the renal excretion of empagliflozin and its glucuronide (EG), the renal glucuronidation of empagliflozin and the urinary glucose excretion.

16-23 % of an oral dose of empagliflozin is excreted unchanged in the urine with a renal clearance of 28-39 ml/min \[Seman2013, Sarashina2013, Macha2014, Macha2014f\]; the glucuronides in the urine are 23 % of the dose \[Chen2015a\].

The renal excretion is relative to the glomerular filtration rate (GFR): the filtration of the unbound fraction (fu 0.14) and a tubular secretion by OAT3 (f\_oat3); probenecid lowers the renal clearance from 32.7 to 15.2 ml/min, the filtration \[Macha2014\].

The kidney glucuronidates empagliflozin (UGT1A9, UGT2B7) after its uptake by OAT3; the renal glucuronidation scales with the renal function (the functional renal mass), which makes the exposure of empagliflozin rise in renal impairment (AUC ratio 1.20, 1.23, 1.67 \[Macha2014f\]).

Empagliflozin lowers the renal threshold of glucose (RTG) by the inhibition of SGLT2; glucose above the RTG is excreted with the GFR.

The content of this model has been carefully created in a manual research effort.

Terms of use

Copyright © 2026 Matthias König.

This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International License.

Redistribution and use of any part of this model, with or without modification, are permitted provided that the following conditions are met:

Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in a different form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Units

| Time | Substance | Extent | Volume | Area | Length |
| --- | --- | --- | --- | --- | --- |
| min | mmol | mmol | l | m<sup>2</sup> | m |

## Compartments

| Symbol | Id | Name | Size | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{Vext}$ | `Vext` | plasma | $1.5$ | l | ✓ |
| $\mathrm{Vki}$ | `Vki` | kidney | $0.3$ | l | ✓ |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{emp}_{\mathrm{ext}}$ | `emp_ext` | empagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{urine}}$ | `emp_urine` | empagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{eg}_{\mathrm{ext}}$ | `eg_ext` | empagliflozin-glucuronide (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{eg}_{\mathrm{urine}}$ | `eg_urine` | empagliflozin-glucuronide (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{emptot}_{\mathrm{urine}}$ | `emptot_urine` | total empagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{fpg}$ | `fpg` | fasting plasma glucose (FPG) | $\mathrm{Vext}$ | $5$ | mmol/l | concentration, boundary |
| $\mathrm{glc}_{\mathrm{urine}}$ | `glc_urine` | glucose (urine) | $\mathrm{Vurine}$ |  | mmol | amount |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $f_{\mathrm{renal\_function}}$ | `f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) | $101$ | ml/min | ✓ |
| $f_{\mathrm{ugt}}$ | `f_ugt` | scaling factor UGT activity | $1$ | \- | ✓ |
| $f_{\mathrm{oat3}}$ | `f_oat3` | scaling factor renal uptake (OAT3) | $1$ | \- | ✓ |
| $\mathrm{fu}_{\mathrm{emp}}$ | `fu_emp` | fraction unbound in plasma empagliflozin | $0.14$ | \- | ✓ |
| $\mathrm{EMPSEC}_{\mathrm{f}}$ | `EMPSEC_f` | tubular secretion of empagliflozin relative to its filtration | $1.31094078052281$ | \- | ✓ |
| $\mathrm{EGEX}_{\mathrm{f}}$ | `EGEX_f` | renal clearance of EG relative to the GFR | $1.66201484620997$ | \- | ✓ |
| $\mathrm{EMP2EG}_{\mathrm{k}}$ | `EMP2EG_k` | rate of the renal glucuronidation of empagliflozin | $0.0566$ | 1/min | ✓ |
| $\mathrm{Mr}_{\mathrm{glc}}$ | `Mr_glc` | molecular weight glucose | $180.16$ | g/mol | ✓ |
| $\mathrm{cf}_{\mathrm{mg\_per\_g}}$ | `cf_mg_per_g` | conversion factor mg per g | $1000$ | mg/g | ✓ |
| $\mathrm{fpg}_{\mathrm{healthy}}$ | `fpg_healthy` | fasting plasma glucose (healthy) | $5$ | mmol/l | ✓ |
| $\mathrm{PG}_{\mathrm{delta}}$ | `PG_delta` | mean plasma glucose of the day minus the fasting plasma glucose | $0.7$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{E50}}$ | `RTG_E50` | EC50 reduction in RTG | $4.08786245088877 \times 10^{-6}$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{gamma}}$ | `RTG_gamma` | hill coefficient reduction in RTG | $1$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{base}}$ | `RTG_base` | baseline RTG | $11.7003562596116$ | mmol/l | ✓ |
| $\mathrm{RTG}_{\mathrm{m\_fpg}}$ | `RTG_m_fpg` | FPG effect on RTG | $2.42163281330897$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{max\_inhibition}}$ | `RTG_max_inhibition` | maximal relative reduction of the RTG | $0.692918965464136$ | \- | ✓ |
| $\mathrm{RTG}_{\mathrm{splay}}$ | `RTG_splay` | splay of the glucose excretion around the RTG | $1$ | mmol/l | ✓ |
| $\mathrm{GFR}$ | `GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{pg}$ | `pg` | mean plasma glucose of the day |  | mmol/l |  |
| $\mathrm{RTG}_{\mathrm{fpg}}$ | `RTG_fpg` | RTG without SGLT2 inhibition |  | mmol/l |  |
| $\mathrm{RTG}_{\mathrm{delta}}$ | `RTG_delta` | maximal reduction of the RTG |  | mmol/l |  |
| $\mathrm{RTG}$ | `RTG` | renal threshold glucose (RTG) |  | mmol/l |  |
| $\mathrm{UGE}$ | `UGE` | urinary glucose excretion (UGE) |  | gram |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{emp}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{eg}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{glc}_{\mathrm{urine}} &= 0 \cdot \mathrm{Vurine}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
\mathrm{GFR} &= f_{\mathrm{renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}} \\[1ex]
\mathrm{emptot}_{\mathrm{urine}} &= \mathrm{emp}_{\mathrm{urine}} + \mathrm{eg}_{\mathrm{urine}} \\[1ex]
\mathrm{pg} &= \mathrm{fpg} + \mathrm{PG}_{\mathrm{delta}} \\[1ex]
\mathrm{RTG}_{\mathrm{fpg}} &= \mathrm{RTG}_{\mathrm{base}} + \mathrm{RTG}_{\mathrm{m\_fpg}} \cdot \mathopen{}\left(\mathrm{fpg} - \mathrm{fpg}_{\mathrm{healthy}}\right) \\[1ex]
\mathrm{RTG}_{\mathrm{delta}} &= \mathrm{RTG}_{\mathrm{fpg}} \cdot \mathrm{RTG}_{\mathrm{max\_inhibition}} \\[1ex]
\mathrm{RTG} &= \mathrm{RTG}_{\mathrm{fpg}} - \frac{\mathrm{RTG}_{\mathrm{delta}} \cdot \mathrm{emp}_{\mathrm{ext}}^{\mathrm{RTG}_{\mathrm{gamma}}}}{\mathrm{RTG}_{\mathrm{E50}}^{\mathrm{RTG}_{\mathrm{gamma}}} + \mathrm{emp}_{\mathrm{ext}}^{\mathrm{RTG}_{\mathrm{gamma}}}} \\[1ex]
\mathrm{UGE} &= \frac{\mathrm{glc}_{\mathrm{urine}} \cdot \mathrm{Mr}_{\mathrm{glc}}}{\mathrm{cf}_{\mathrm{mg\_per\_g}}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{EMP2EG}}$ | `EMP2EG` | empagliflozin glucuronidation (EMP2EG) | $\mathrm{emp}_{\mathrm{ext}} \longrightarrow \mathrm{eg}_{\mathrm{ext}}$ |  |
| $v_{\mathrm{EMPEX}}$ | `EMPEX` | empagliflozin excretion (EMPEX) | $\mathrm{emp}_{\mathrm{ext}} \longrightarrow \mathrm{emp}_{\mathrm{urine}}$ |  |
| $v_{\mathrm{EGEX}}$ | `EGEX` | EG excretion (EGEX) | $\mathrm{eg}_{\mathrm{ext}} \longrightarrow \mathrm{eg}_{\mathrm{urine}}$ |  |
| $v_{\mathrm{GLCEX}}$ | `GLCEX` | glucose excretion (GLCEX) | $\mathrm{fpg} \longrightarrow \mathrm{glc}_{\mathrm{urine}}$ | $\mathrm{emp}_{\mathrm{ext}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{EMP2EG}} &= f_{\mathrm{renal\_function}} \cdot f_{\mathrm{oat3}} \cdot f_{\mathrm{ugt}} \cdot \mathrm{EMP2EG}_{\mathrm{k}} \cdot \mathrm{Vki} \cdot \mathrm{emp}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{EMPEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{fu}_{\mathrm{emp}} \cdot \mathopen{}\left(1 + f_{\mathrm{oat3}} \cdot \mathrm{EMPSEC}_{\mathrm{f}}\right) \cdot \mathrm{emp}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{EGEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{EGEX}_{\mathrm{f}} \cdot \mathrm{eg}_{\mathrm{ext}} \\[1ex]
v_{\mathrm{GLCEX}} &= \mathrm{GFR} \cdot 0.001 \cdot \mathrm{RTG}_{\mathrm{splay}} \cdot \ln\mathopen{}\left(1 + \exp\mathopen{}\left(\frac{\mathrm{pg} - \mathrm{RTG}}{\mathrm{RTG}_{\mathrm{splay}}}\right)\right)
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{EMP2EG}} - v_{\mathrm{EMPEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{EMPEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{eg}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMP2EG}} - v_{\mathrm{EGEX}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{eg}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{EGEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{glc}_{\mathrm{urine}}}{\mathrm{d} t} &= v_{\mathrm{GLCEX}}
\end{aligned}
$$
