# Model for the hepatic metabolism of empagliflozin.

Model `empagliflozin_liver`, SBML Level 3 Version 2, read from empagliflozin\_liver.xml, written by sbmlode 0.3.0.

Model for the hepatic metabolism of empagliflozin.

Empagliflozin is taken up into the hepatocytes by OATP1B1 and OATP1B3 and glucuronidated by UGT2B7, UGT1A3, UGT1A8 and UGT1A9 to the 2-O-, 3-O- and 6-O-glucuronides (EG), which are released into the plasma and excreted in the urine; no glucuronide is found in the feces \[Chen2015a\].

The hepatic clearance is a first-order clearance of the plasma of the liver: the uptake is the rate determining step of the low permeable empagliflozin (f\_oatp), the glucuronidation (f\_ugt) is far below its Km. It scales with the volume of the liver tissue, so with the tissue loss in cirrhosis.

The exposure rises with the Child-Pugh class (AUC ratio 1.23, 1.47, 1.75 \[Macha2014b\]); with rifampicin (single dose, inhibition of OATP1B1/1B3) the AUC rises by 35 % \[Macha2014\].

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
| $\mathrm{Vli}$ | `Vli` | liver | $1.5$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{emp}_{\mathrm{ext}}$ | `emp_ext` | empagliflozin (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |
| $\mathrm{eg}_{\mathrm{ext}}$ | `eg_ext` | empagliflozin-glucuronide (plasma) | $\mathrm{Vext}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{EMP2EG}_{\mathrm{k}}$ | `EMP2EG_k` | rate of the hepatic glucuronidation of empagliflozin | $0.0176461728190058$ | 1/min | ✓ |
| $f_{\mathrm{ugt}}$ | `f_ugt` | scaling factor UGT activity | $1$ | \- | ✓ |
| $f_{\mathrm{oatp}}$ | `f_oatp` | scaling factor hepatic uptake (OATP1B1/1B3) | $1$ | \- | ✓ |

## Reactions

| Rate | Id | Name | Equation |
| --- | --- | --- | --- |
| $v_{\mathrm{EMP2EG}}$ | `EMP2EG` | empagliflozin glucuronidation (EMP2EG) | $\mathrm{emp}_{\mathrm{ext}} \longrightarrow \mathrm{eg}_{\mathrm{ext}}$ |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{EMP2EG}} &= f_{\mathrm{oatp}} \cdot f_{\mathrm{ugt}} \cdot \mathrm{EMP2EG}_{\mathrm{k}} \cdot \mathrm{Vli} \cdot \mathrm{emp}_{\mathrm{ext}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{emp}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{EMP2EG}}}{\mathrm{Vext}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{eg}_{\mathrm{ext}}}{\mathrm{d} t} &= \frac{v_{\mathrm{EMP2EG}}}{\mathrm{Vext}}
\end{aligned}
$$
