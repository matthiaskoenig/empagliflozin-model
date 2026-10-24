# Model for empagliflozin absorption in the small intestine

Model `empagliflozin_intestine`, SBML Level 3 Version 2, read from empagliflozin\_intestine.xml, written by sbmlode 0.3.0.

Model for empagliflozin absorption

Dissolution of the oral dose and absorption of empagliflozin; the fraction absorbed F\_emp\_abs of the empagliflozin of the lumen is absorbed, the rest is excreted in the feces after the transit through the intestine (METEXC\_k). After an oral solution of 50 mg of 14C-empagliflozin 41.2 % of the radioactivity is recovered in the feces, 34.2 % of the dose as unchanged empagliflozin, no glucuronide; 54.4 % in the urine \[Chen2015a\]. The fecal excretion is slow: 10.9, 29.0, 37.9 and 39.9 % of the dose at 48, 72, 96 and 120 hr.

Food (fed) slows the absorption (FED\_f\_absorption): a high fat breakfast lowers Cmax by 29 % and delays tmax from 1.0 to 2.5 hr, the AUC is lowered by 10 % \[Seman2013\]; a Chinese high fat meal lowers Cmax by 20 % and the AUC by 5 % \[Chen2020\].

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
| $\mathrm{Vext}$ | `Vext` | plasma | $1$ | l | ✓ |
| $\mathrm{Vgu}$ | `Vgu` | intestine | $1.2825$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) | $1.15425$ | l |  |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vapical}$ | `Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |
| $\mathrm{Vintestine}_{0}$ | `Vintestine_0` | intestinal lumen (transit 0) | $0.23085$ | l | ✓ |
| $\mathrm{Vintestine}_{1}$ | `Vintestine_1` | intestinal lumen (transit 1) | $0.23085$ | l | ✓ |
| $\mathrm{Vintestine}_{2}$ | `Vintestine_2` | intestinal lumen (transit 2) | $0.23085$ | l | ✓ |
| $\mathrm{Vintestine}_{3}$ | `Vintestine_3` | intestinal lumen (transit 3) | $0.23085$ | l | ✓ |
| $\mathrm{Vintestine}_{4}$ | `Vintestine_4` | intestinal lumen (transit 4) | $0.23085$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{emp}_{\mathrm{stomach}}$ | `emp_stomach` | empagliflozin (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{emp}_{\mathrm{lumen}}$ | `emp_lumen` | empagliflozin (intestinal lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{ext}}$ | `emp_ext` | empagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{feces}}$ | `emp_feces` | empagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{emp}_{\mathrm{intestine\_0}}$ | `emp_intestine_0` | empagliflozin (intestinal transit 0) | $\mathrm{Vintestine}_{0}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{intestine\_1}}$ | `emp_intestine_1` | empagliflozin (intestinal transit 1) | $\mathrm{Vintestine}_{1}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{intestine\_2}}$ | `emp_intestine_2` | empagliflozin (intestinal transit 2) | $\mathrm{Vintestine}_{2}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{intestine\_3}}$ | `emp_intestine_3` | empagliflozin (intestinal transit 3) | $\mathrm{Vintestine}_{3}$ | $0$ | mmol/l | concentration |
| $\mathrm{emp}_{\mathrm{intestine\_4}}$ | `emp_intestine_4` | empagliflozin (intestinal transit 4) | $\mathrm{Vintestine}_{4}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $F_{\mathrm{emp\_abs}}$ | `F_emp_abs` | fraction absorbed empagliflozin | $0.438440610233251$ | \- | ✓ |
| $\mathrm{EMPABS}_{\mathrm{k}}$ | `EMPABS_k` | rate of empagliflozin absorption | $0.0129016148159033$ | 1/min | ✓ |
| $\mathrm{METEXC}_{\mathrm{k}}$ | `METEXC_k` | rate of the fecal transit | $0.000785814929207019$ | 1/min | ✓ |
| $\mathrm{fed}$ | `fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{FED}_{\mathrm{f\_absorption}}$ | `FED_f_absorption` | fed effect on the absorption rate | $0.549647454138511$ | \- | ✓ |
| $f_{\mathrm{absorption}}$ | `f_absorption` | factor of the absorption rate |  | \- |  |
| $\mathrm{PODOSE}_{\mathrm{emp}}$ | `PODOSE_emp` | oral dose empagliflozin \[mg\] | $0$ | mg |  |
| $\mathrm{Ka}_{\mathrm{dis\_emp}}$ | `Ka_dis_emp` | dissolution rate empagliflozin | $2$ | 1/hr | ✓ |
| $\mathrm{Mr}_{\mathrm{emp}}$ | `Mr_emp` | Molecular weight empagliflozin \[g/mole\] | $450.909$ | g/mol | ✓ |
| $\mathrm{absorption}$ | `absorption` | empagliflozin leaving the lumen |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{emp}_{\mathrm{stomach}} &= 0 \cdot \mathrm{Vstomach} \\[1ex]
\mathrm{emp}_{\mathrm{feces}} &= 0 \cdot \mathrm{Vfeces}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
f_{\mathrm{absorption}} &= 1 + \mathrm{fed} \cdot \mathopen{}\left(\mathrm{FED}_{\mathrm{f\_absorption}} - 1\right) \\[1ex]
\mathrm{absorption} &= f_{\mathrm{absorption}} \cdot \mathrm{EMPABS}_{\mathrm{k}} \cdot \mathrm{Vlumen} \cdot \mathrm{emp}_{\mathrm{lumen}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{EMPABS}}$ | `EMPABS` | absorption empagliflozin (EMPABS) | $\mathrm{emp}_{\mathrm{lumen}} \longrightarrow \mathrm{emp}_{\mathrm{ext}}$ |
| $v_{\mathrm{EMPEXC}}$ | `EMPEXC` | empagliflozin to the fecal transit (EMPEXC) | $\mathrm{emp}_{\mathrm{lumen}} \longrightarrow \mathrm{emp}_{\mathrm{intestine\_0}}$ |
| $v_{\mathrm{EMPEXC\_0}}$ | `EMPEXC_0` | fecal transit empagliflozin 0 | $\mathrm{emp}_{\mathrm{intestine\_0}} \longrightarrow \mathrm{emp}_{\mathrm{intestine\_1}}$ |
| $v_{\mathrm{EMPEXC\_1}}$ | `EMPEXC_1` | fecal transit empagliflozin 1 | $\mathrm{emp}_{\mathrm{intestine\_1}} \longrightarrow \mathrm{emp}_{\mathrm{intestine\_2}}$ |
| $v_{\mathrm{EMPEXC\_2}}$ | `EMPEXC_2` | fecal transit empagliflozin 2 | $\mathrm{emp}_{\mathrm{intestine\_2}} \longrightarrow \mathrm{emp}_{\mathrm{intestine\_3}}$ |
| $v_{\mathrm{EMPEXC\_3}}$ | `EMPEXC_3` | fecal transit empagliflozin 3 | $\mathrm{emp}_{\mathrm{intestine\_3}} \longrightarrow \mathrm{emp}_{\mathrm{intestine\_4}}$ |
| $v_{\mathrm{EMPEXC\_4}}$ | `EMPEXC_4` | fecal transit empagliflozin 4 | $\mathrm{emp}_{\mathrm{intestine\_4}} \longrightarrow \mathrm{emp}_{\mathrm{feces}}$ |
| $v_{\mathrm{dissolution\_emp}}$ | `dissolution_emp` | dissolution empagliflozin | $\mathrm{emp}_{\mathrm{stomach}} \longrightarrow \mathrm{emp}_{\mathrm{lumen}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{EMPABS}} &= F_{\mathrm{emp\_abs}} \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{EMPEXC}} &= \mathopen{}\left(1 - F_{\mathrm{emp\_abs}}\right) \cdot \mathrm{absorption} \\[1ex]
v_{\mathrm{EMPEXC\_0}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vintestine}_{0} \cdot \mathrm{emp}_{\mathrm{intestine\_0}} \\[1ex]
v_{\mathrm{EMPEXC\_1}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vintestine}_{1} \cdot \mathrm{emp}_{\mathrm{intestine\_1}} \\[1ex]
v_{\mathrm{EMPEXC\_2}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vintestine}_{2} \cdot \mathrm{emp}_{\mathrm{intestine\_2}} \\[1ex]
v_{\mathrm{EMPEXC\_3}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vintestine}_{3} \cdot \mathrm{emp}_{\mathrm{intestine\_3}} \\[1ex]
v_{\mathrm{EMPEXC\_4}} &= \mathrm{METEXC}_{\mathrm{k}} \cdot \mathrm{Vintestine}_{4} \cdot \mathrm{emp}_{\mathrm{intestine\_4}} \\[1ex]
v_{\mathrm{dissolution\_emp}} &= \frac{\frac{\mathrm{Ka}_{\mathrm{dis\_emp}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{emp}}}{\mathrm{Mr}_{\mathrm{emp}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{lumen}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{EMPABS}} - v_{\mathrm{EMPEXC}} + v_{\mathrm{dissolution\_emp}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPABS}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{feces}}}{\mathrm{d} t} &= v_{\mathrm{EMPEXC\_4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{intestine\_0}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPEXC}} - v_{\mathrm{EMPEXC\_0}}}{\mathrm{Vintestine}_{0}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{intestine\_1}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPEXC\_0}} - v_{\mathrm{EMPEXC\_1}}}{\mathrm{Vintestine}_{1}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{intestine\_2}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPEXC\_1}} - v_{\mathrm{EMPEXC\_2}}}{\mathrm{Vintestine}_{2}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{intestine\_3}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPEXC\_2}} - v_{\mathrm{EMPEXC\_3}}}{\mathrm{Vintestine}_{3}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{intestine\_4}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMPEXC\_3}} - v_{\mathrm{EMPEXC\_4}}}{\mathrm{Vintestine}_{4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{emp}}}{\mathrm{d} t} &= -v_{\mathrm{dissolution\_emp}} \cdot \mathrm{Mr}_{\mathrm{emp}} \qquad \text{(rate rule)}
\end{aligned}
$$
