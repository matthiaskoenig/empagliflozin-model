# Empagliflozin body model

Model `empagliflozin_body`, SBML Level 3 Version 2, read from empagliflozin\_body.xml, written by sbmlode 0.3.0.

Empagliflozin body model

Distribution of empagliflozin and its glucuronide (EG) with the blood flow, the tissue distribution of empagliflozin into a fast and a deep pool, the hepatic and renal glucuronidation, the absorption and fecal excretion and the renal excretion of empagliflozin and EG.

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
| $\mathrm{Vre}$ | `Vre` | rest of body |  | l |  |
| $\mathrm{Vgu}$ | `Vgu` | gut |  | l |  |
| $\mathrm{Vki}$ | `Vki` | kidney |  | l |  |
| $\mathrm{Vli}$ | `Vli` | liver |  | l |  |
| $\mathrm{Vlu}$ | `Vlu` | lung |  | l |  |
| $\mathrm{Vve}$ | `Vve` | venous blood |  | l |  |
| $\mathrm{Var}$ | `Var` | arterial blood |  | l |  |
| $\mathrm{Vurine}$ | `Vurine` | urine | $1$ | l | ✓ |
| $\mathrm{Vfeces}$ | `Vfeces` | feces | $1$ | l | ✓ |
| $\mathrm{Vstomach}$ | `Vstomach` | stomach | $1$ | l | ✓ |
| $\mathrm{Vlumen}$ | `Vlumen` | intestinal lumen (inner part of intestine) |  | l |  |
| $\mathrm{Vpo}$ | `Vpo` | portal plasma |  | l |  |
| $\mathrm{Vhv}$ | `Vhv` | hepatic venous plasma |  | l |  |
| $\mathrm{Vki}_{\mathrm{tissue}}$ | `Vki_tissue` | ki tissue |  | l |  |
| $\mathrm{Vki}_{\mathrm{plasma}}$ | `Vki_plasma` | ki plasma |  | l |  |
| $\mathrm{Vli}_{\mathrm{tissue}}$ | `Vli_tissue` | li tissue |  | l |  |
| $\mathrm{Vli}_{\mathrm{plasma}}$ | `Vli_plasma` | li plasma |  | l |  |
| $\mathrm{Vlu}_{\mathrm{tissue}}$ | `Vlu_tissue` | lu tissue |  | l |  |
| $\mathrm{Vlu}_{\mathrm{plasma}}$ | `Vlu_plasma` | lu plasma |  | l |  |
| $\mathrm{Vgu}_{\mathrm{tissue}}$ | `Vgu_tissue` | gu tissue |  | l |  |
| $\mathrm{Vgu}_{\mathrm{plasma}}$ | `Vgu_plasma` | gu plasma |  | l |  |
| $\mathrm{Vre}_{\mathrm{tissue}}$ | `Vre_tissue` | re tissue |  | l |  |
| $\mathrm{Vre}_{\mathrm{plasma}}$ | `Vre_plasma` | re plasma |  | l |  |
| $\mathrm{GU}_{\mathrm{\_Vapical}}$ | `GU__Vapical` | apical membrane (intestinal membrane enterocytes) | $\mathrm{NaN}$ | m<sup>2</sup> | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vintestine\_0}}$ | `GU__Vintestine_0` | intestinal lumen (transit 0) | $0.23085$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vintestine\_1}}$ | `GU__Vintestine_1` | intestinal lumen (transit 1) | $0.23085$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vintestine\_2}}$ | `GU__Vintestine_2` | intestinal lumen (transit 2) | $0.23085$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vintestine\_3}}$ | `GU__Vintestine_3` | intestinal lumen (transit 3) | $0.23085$ | l | ✓ |
| $\mathrm{GU}_{\mathrm{\_Vintestine\_4}}$ | `GU__Vintestine_4` | intestinal lumen (transit 4) | $0.23085$ | l | ✓ |

## Species

| Symbol | Id | Name | Compartment | Value | Unit | Properties |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{Cki}_{\mathrm{plasma\_emp}}$ | `Cki_plasma_emp` | empagliflozin (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_emp}}$ | `Cli_plasma_emp` | empagliflozin (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_emp}}$ | `Clu_plasma_emp` | empagliflozin (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_emp}}$ | `Cgu_plasma_emp` | empagliflozin (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_emp}}$ | `Cre_plasma_emp` | empagliflozin (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{emp}}$ | `Car_emp` | empagliflozin (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{emp}}$ | `Cve_emp` | empagliflozin (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{emp}}$ | `Cpo_emp` | empagliflozin (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{emp}}$ | `Chv_emp` | empagliflozin (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{emp}}$ | `Clu_emp` | empagliflozin (lung) | $\mathrm{Vlu}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{emp}}$ | `Cre_emp` | empagliflozin (rest) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{deep\_emp}}$ | `Cre_deep_emp` | empagliflozin (rest, deep tissue) | $\mathrm{Vre}_{\mathrm{tissue}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{emp}}$ | `Aurine_emp` | empagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Cki}_{\mathrm{plasma\_eg}}$ | `Cki_plasma_eg` | empagliflozin-glucuronide (kidney plasma) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cli}_{\mathrm{plasma\_eg}}$ | `Cli_plasma_eg` | empagliflozin-glucuronide (liver plasma) | $\mathrm{Vli}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Clu}_{\mathrm{plasma\_eg}}$ | `Clu_plasma_eg` | empagliflozin-glucuronide (lung plasma) | $\mathrm{Vlu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cgu}_{\mathrm{plasma\_eg}}$ | `Cgu_plasma_eg` | empagliflozin-glucuronide (gut plasma) | $\mathrm{Vgu}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cre}_{\mathrm{plasma\_eg}}$ | `Cre_plasma_eg` | empagliflozin-glucuronide (rest plasma) | $\mathrm{Vre}_{\mathrm{plasma}}$ | $0$ | mmol/l | concentration |
| $\mathrm{Car}_{\mathrm{eg}}$ | `Car_eg` | empagliflozin-glucuronide (arterial blood plasma) | $\mathrm{Var}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cve}_{\mathrm{eg}}$ | `Cve_eg` | empagliflozin-glucuronide (venous blood plasma) | $\mathrm{Vve}$ | $0$ | mmol/l | concentration |
| $\mathrm{Cpo}_{\mathrm{eg}}$ | `Cpo_eg` | empagliflozin-glucuronide (portal vein plasma) | $\mathrm{Vpo}$ | $0$ | mmol/l | concentration |
| $\mathrm{Chv}_{\mathrm{eg}}$ | `Chv_eg` | empagliflozin-glucuronide (hepatic vein plasma) | $\mathrm{Vhv}$ | $0$ | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{eg}}$ | `Aurine_eg` | empagliflozin-glucuronide (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Cgu}_{\mathrm{emp}}$ | `Cgu_emp` | empagliflozin (intestinal lumen) | $\mathrm{Vlumen}$ | $0$ | mmol/l | concentration |
| $\mathrm{Afeces}_{\mathrm{emp}}$ | `Afeces_emp` | empagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{Cve}_{\mathrm{emptot}}$ | `Cve_emptot` | total empagliflozin (plasma) | $\mathrm{Vve}$ |  | mmol/l | concentration |
| $\mathrm{Aurine}_{\mathrm{emptot}}$ | `Aurine_emptot` | total empagliflozin (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{Afeces}_{\mathrm{emptot}}$ | `Afeces_emptot` | total empagliflozin (feces) | $\mathrm{Vfeces}$ |  | mmol | amount |
| $\mathrm{KI}_{\mathrm{\_fpg}}$ | `KI__fpg` | fasting plasma glucose (FPG) | $\mathrm{Vki}_{\mathrm{plasma}}$ | $5$ | mmol/l | concentration, boundary |
| $\mathrm{KI}_{\mathrm{\_glc\_urine}}$ | `KI__glc_urine` | glucose (urine) | $\mathrm{Vurine}$ |  | mmol | amount |
| $\mathrm{GU}_{\mathrm{\_emp\_stomach}}$ | `GU__emp_stomach` | empagliflozin (stomach) | $\mathrm{Vstomach}$ |  | mmol | amount, boundary |
| $\mathrm{GU}_{\mathrm{\_emp\_intestine\_0}}$ | `GU__emp_intestine_0` | empagliflozin (intestinal transit 0) | $\mathrm{GU}_{\mathrm{\_Vintestine\_0}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_emp\_intestine\_1}}$ | `GU__emp_intestine_1` | empagliflozin (intestinal transit 1) | $\mathrm{GU}_{\mathrm{\_Vintestine\_1}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_emp\_intestine\_2}}$ | `GU__emp_intestine_2` | empagliflozin (intestinal transit 2) | $\mathrm{GU}_{\mathrm{\_Vintestine\_2}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_emp\_intestine\_3}}$ | `GU__emp_intestine_3` | empagliflozin (intestinal transit 3) | $\mathrm{GU}_{\mathrm{\_Vintestine\_3}}$ | $0$ | mmol/l | concentration |
| $\mathrm{GU}_{\mathrm{\_emp\_intestine\_4}}$ | `GU__emp_intestine_4` | empagliflozin (intestinal transit 4) | $\mathrm{GU}_{\mathrm{\_Vintestine\_4}}$ | $0$ | mmol/l | concentration |

## Parameters

| Symbol | Id | Name | Value | Unit | Constant |
| --- | --- | --- | --- | --- | :---: |
| $\mathrm{BW}$ | `BW` | body weight | $75$ | kg | ✓ |
| $\mathrm{BW}_{\mathrm{ref}}$ | `BW_ref` | reference body weight \[kg\] | $75$ | kg | ✓ |
| $\mathrm{COBW}$ | `COBW` | cardiac output per bodyweight | $1.548$ | ml/s/kg | ✓ |
| $f_{\mathrm{cardiac\_function}}$ | `f_cardiac_function` | heart function | $1$ | \- |  |
| $\mathrm{CO}$ | `CO` | cardiac output \[ml/s\] |  | ml/s |  |
| $\mathrm{QC}$ | `QC` | cardiac output \[l/min\] |  | l/min |  |
| $\mathrm{Fblood}$ | `Fblood` | blood fraction of organ volume | $0.02$ | \- |  |
| $\mathrm{HCT}$ | `HCT` | hematocrit | $0.51$ | \- | ✓ |
| $\mathrm{FVgu}$ | `FVgu` | gut fractional tissue volume | $0.0171$ | l/kg | ✓ |
| $\mathrm{FVki}$ | `FVki` | kidney fractional tissue volume | $0.0044$ | l/kg | ✓ |
| $\mathrm{FVli}$ | `FVli` | liver fractional tissue volume | $0.021$ | l/kg | ✓ |
| $\mathrm{FVlu}$ | `FVlu` | lung fractional tissue volume | $0.0076$ | l/kg | ✓ |
| $\mathrm{FVre}$ | `FVre` | rest of body fractional tissue volume |  | l/kg |  |
| $\mathrm{FVve}$ | `FVve` | venous fractional tissue volume | $0.0514$ | l/kg | ✓ |
| $\mathrm{FVar}$ | `FVar` | arterial fractional tissue volume | $0.0257$ | l/kg | ✓ |
| $\mathrm{FVpo}$ | `FVpo` | portal fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FVhv}$ | `FVhv` | hepatic venous fractional tissue volume | $0.001$ | l/kg | ✓ |
| $\mathrm{FQgu}$ | `FQgu` | gut fractional tissue blood flow | $0.18$ | \- | ✓ |
| $\mathrm{FQki}$ | `FQki` | kidney fractional tissue blood flow | $0.19$ | \- | ✓ |
| $\mathrm{FQh}$ | `FQh` | hepatic (venous side) fractional tissue blood flow | $0.215$ | \- | ✓ |
| $\mathrm{FQlu}$ | `FQlu` | lung fractional tissue blood flow | $1$ | \- | ✓ |
| $\mathrm{FQre}$ | `FQre` | rest of body fractional tissue blood flow |  | \- |  |
| $f_{\mathrm{cirrhosis}}$ | `f_cirrhosis` | severity of cirrhosis \[0, 0.95\] | $0$ | \- | ✓ |
| $f_{\mathrm{shunts}}$ | `f_shunts` | fraction of portal venous blood shunted by the liver |  | \- |  |
| $f_{\mathrm{tissue\_loss}}$ | `f_tissue_loss` | fraction of lost parenchymal liver volume |  | \- |  |
| $f_{\mathrm{GFR\_Qki}}$ | `f_GFR_Qki` | relationship GFR and Qki | $76.3105006271061$ | ml/l | ✓ |
| $\mathrm{GFR}_{\mathrm{base}}$ | `GFR_base` | glomerular filtration rate (healthy) |  | ml/min |  |
| $f_{\mathrm{ugt}}$ | `f_ugt` | scaling factor UGT activity | $1$ | \- | ✓ |
| $\mathrm{EMP2EG}_{\mathrm{k}}$ | `EMP2EG_k` | rate of the hepatic glucuronidation of empagliflozin | $0.0176461728190058$ | 1/min | ✓ |
| $f_{\mathrm{ugt\_renal}}$ | `f_ugt_renal` | renal fraction of the glucuronidation | $0.3$ | \- | ✓ |
| $\mathrm{EMP2EG}_{\mathrm{k\_renal}}$ | `EMP2EG_k_renal` | rate of the renal glucuronidation of empagliflozin |  | 1/min |  |
| $\mathrm{PODOSE}_{\mathrm{emp}}$ | `PODOSE_emp` | oral dose emp \[mg\] | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{emp}}$ | `Mr_emp` | Molecular weight emp \[g/mole\] | $450.909$ | g/mol | ✓ |
| $\mathrm{ftissue}_{\mathrm{emp}}$ | `ftissue_emp` | tissue distribution emp | $1.11016781396837$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{emp}}$ | `Kp_emp` | tissue/plasma partition coefficient emp | $0.238160000879157$ | \- | ✓ |
| $\mathrm{ftissue}_{\mathrm{deep\_emp}}$ | `ftissue_deep_emp` | deep tissue distribution emp | $0.112608734639073$ | l/min | ✓ |
| $\mathrm{Kp}_{\mathrm{deep\_emp}}$ | `Kp_deep_emp` | deep tissue/plasma partition coefficient emp | $0.418946866479524$ | \- | ✓ |
| $\mathrm{IVDOSE}_{\mathrm{emp}}$ | `IVDOSE_emp` | IV bolus dose emp \[mg\] | $0$ | mg |  |
| $\mathrm{ti}_{\mathrm{emp}}$ | `ti_emp` | injection time emp \[s\] | $10$ | s | ✓ |
| $\mathrm{Ki}_{\mathrm{emp}}$ | `Ki_emp` | Ki \[1/min\] injection emp |  | 1/min |  |
| $\mathrm{Ri}_{\mathrm{emp}}$ | `Ri_emp` | Ri \[mg/min\] rate of infusion emp | $0$ | mg/min | ✓ |
| $\mathrm{cum}_{\mathrm{dose\_emp}}$ | `cum_dose_emp` | Cumulative dose due to infusion emp | $0$ | mg |  |
| $\mathrm{Mr}_{\mathrm{eg}}$ | `Mr_eg` | Molecular weight eg \[g/mole\] | $627$ | g/mol | ✓ |
| $\mathrm{Qgu}$ | `Qgu` | gut blood flow |  | l/min |  |
| $\mathrm{Qki}$ | `Qki` | kidney blood flow |  | l/min |  |
| $\mathrm{Qh}$ | `Qh` | hepatic (venous side) blood flow |  | l/min |  |
| $\mathrm{Qha}$ | `Qha` | hepatic artery blood flow |  | l/min |  |
| $\mathrm{Qlu}$ | `Qlu` | lung blood flow |  | l/min |  |
| $\mathrm{Qre}$ | `Qre` | rest of body blood flow |  | l/min |  |
| $\mathrm{Qpo}$ | `Qpo` | portal blood flow |  | l/min |  |
| $\mathrm{KI}_{\mathrm{\_f\_renal\_function}}$ | `KI__f_renal_function` | renal function | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_f\_oat3}}$ | `KI__f_oat3` | scaling factor renal uptake (OAT3) | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_fu\_emp}}$ | `KI__fu_emp` | fraction unbound in plasma empagliflozin | $0.14$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_EMPSEC\_f}}$ | `KI__EMPSEC_f` | tubular secretion of empagliflozin relative to its filtration | $1.31094078052281$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_EGEX\_f}}$ | `KI__EGEX_f` | renal clearance of EG relative to the GFR | $1.66201484620997$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_Mr\_glc}}$ | `KI__Mr_glc` | molecular weight glucose | $180.16$ | g/mol | ✓ |
| $\mathrm{KI}_{\mathrm{\_cf\_mg\_per\_g}}$ | `KI__cf_mg_per_g` | conversion factor mg per g | $1000$ | mg/g | ✓ |
| $\mathrm{KI}_{\mathrm{\_fpg\_healthy}}$ | `KI__fpg_healthy` | fasting plasma glucose (healthy) | $5$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_PG\_delta}}$ | `KI__PG_delta` | mean plasma glucose of the day minus the fasting plasma glucose | $0.7$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_E50}}$ | `KI__RTG_E50` | EC50 reduction in RTG | $4.08786245088877 \times 10^{-6}$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_gamma}}$ | `KI__RTG_gamma` | hill coefficient reduction in RTG | $1$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_base}}$ | `KI__RTG_base` | baseline RTG | $11.7003562596116$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_m\_fpg}}$ | `KI__RTG_m_fpg` | FPG effect on RTG | $2.42163281330897$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_max\_inhibition}}$ | `KI__RTG_max_inhibition` | maximal relative reduction of the RTG | $0.692918965464136$ | \- | ✓ |
| $\mathrm{KI}_{\mathrm{\_RTG\_splay}}$ | `KI__RTG_splay` | splay of the glucose excretion around the RTG | $1$ | mmol/l | ✓ |
| $\mathrm{KI}_{\mathrm{\_GFR}}$ | `KI__GFR` | glomerular filtration rate |  | ml/min |  |
| $\mathrm{KI}_{\mathrm{\_pg}}$ | `KI__pg` | mean plasma glucose of the day |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG\_fpg}}$ | `KI__RTG_fpg` | RTG without SGLT2 inhibition |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG\_delta}}$ | `KI__RTG_delta` | maximal reduction of the RTG |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_RTG}}$ | `KI__RTG` | renal threshold glucose (RTG) |  | mmol/l |  |
| $\mathrm{KI}_{\mathrm{\_UGE}}$ | `KI__UGE` | urinary glucose excretion (UGE) |  | gram |  |
| $\mathrm{LI}_{\mathrm{\_f\_oatp}}$ | `LI__f_oatp` | scaling factor hepatic uptake (OATP1B1/1B3) | $1$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_F\_emp\_abs}}$ | `GU__F_emp_abs` | fraction absorbed empagliflozin | $0.438440610233251$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_EMPABS\_k}}$ | `GU__EMPABS_k` | rate of empagliflozin absorption | $0.0129016148159033$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_METEXC\_k}}$ | `GU__METEXC_k` | rate of the fecal transit | $0.000785814929207019$ | 1/min | ✓ |
| $\mathrm{GU}_{\mathrm{\_fed}}$ | `GU__fed` | fed state (0 fasted, 1 fed) | $0$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}}$ | `GU__FED_f_absorption` | fed effect on the absorption rate | $0.549647454138511$ | \- | ✓ |
| $\mathrm{GU}_{\mathrm{\_f\_absorption}}$ | `GU__f_absorption` | factor of the absorption rate |  | \- |  |
| $\mathrm{GU}_{\mathrm{\_Ka\_dis\_emp}}$ | `GU__Ka_dis_emp` | dissolution rate empagliflozin | $2$ | 1/hr | ✓ |
| $\mathrm{GU}_{\mathrm{\_absorption}}$ | `GU__absorption` | empagliflozin leaving the lumen |  | mmol/min |  |

## Initial assignments and assignment rules

The initial assignments set the values at $t = 0$:

$$
\begin{aligned}
\mathrm{Aurine}_{\mathrm{emp}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Aurine}_{\mathrm{eg}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{Afeces}_{\mathrm{emp}} &= 0 \cdot \mathrm{Vfeces} \\[1ex]
\mathrm{KI}_{\mathrm{\_glc\_urine}} &= 0 \cdot \mathrm{Vurine} \\[1ex]
\mathrm{GU}_{\mathrm{\_emp\_stomach}} &= 0 \cdot \mathrm{Vstomach}
\end{aligned}
$$

The assignment rules hold at every time $t$:

$$
\begin{aligned}
f_{\mathrm{shunts}} &= f_{\mathrm{cirrhosis}} \\[1ex]
f_{\mathrm{tissue\_loss}} &= f_{\mathrm{cirrhosis}} \\[1ex]
\mathrm{Cve}_{\mathrm{emptot}} &= \mathrm{Cve}_{\mathrm{emp}} + \mathrm{Cve}_{\mathrm{eg}} \\[1ex]
\mathrm{Afeces}_{\mathrm{emptot}} &= \mathrm{Afeces}_{\mathrm{emp}} \\[1ex]
\mathrm{EMP2EG}_{\mathrm{k\_renal}} &= \frac{\frac{\mathrm{EMP2EG}_{\mathrm{k}} \cdot \mathrm{FVli}}{\mathrm{FVki}} \cdot f_{\mathrm{ugt\_renal}}}{1 - f_{\mathrm{ugt\_renal}}} \\[1ex]
\mathrm{FVre} &= 1 - \mathrm{FVgu} - \mathrm{FVki} - \mathrm{FVli} \\
&\quad {} - \mathrm{FVlu} - \mathrm{FVve} - \mathrm{FVar} \\[1ex]
\mathrm{FQre} &= 1 - \mathrm{FQki} - \mathrm{FQh} \\[1ex]
\mathrm{CO} &= f_{\mathrm{cardiac\_function}} \cdot \mathrm{BW} \cdot \mathrm{COBW} \\[1ex]
\mathrm{QC} &= \frac{\mathrm{CO}}{1000} \cdot 60 \\[1ex]
\mathrm{Vgu} &= \mathrm{BW} \cdot \mathrm{FVgu} \\[1ex]
\mathrm{Vki} &= \mathrm{BW} \cdot \mathrm{FVki} \\[1ex]
\mathrm{Vli} &= \mathrm{BW} \cdot \mathrm{FVli} \\[1ex]
\mathrm{Vlu} &= \mathrm{BW} \cdot \mathrm{FVlu} \\[1ex]
\mathrm{Vre} &= \mathrm{BW} \cdot \mathrm{FVre} \\[1ex]
\mathrm{Vve} &= \mathrm{BW} \cdot \mathrm{FVve} - \frac{\mathrm{FVve}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Var} &= \mathrm{BW} \cdot \mathrm{FVar} - \frac{\mathrm{FVar}}{\mathrm{FVar} + \mathrm{FVve}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVve} - \mathrm{FVar}\right) \\[1ex]
\mathrm{Vpo} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVpo} - \frac{\mathrm{FVpo}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Vhv} &= \mathopen{}\left(1 - \mathrm{HCT}\right) \cdot \mathopen{}\left(\mathrm{BW} \cdot \mathrm{FVhv} - \frac{\mathrm{FVhv}}{\mathrm{FVar} + \mathrm{FVve} + \mathrm{FVpo} + \mathrm{FVhv}} \cdot \mathrm{BW} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{FVar} - \mathrm{FVve} - \mathrm{FVpo} - \mathrm{FVhv}\right)\right) \\[1ex]
\mathrm{Qgu} &= \mathrm{QC} \cdot \mathrm{FQgu} \\[1ex]
\mathrm{Qki} &= \mathrm{QC} \cdot \mathrm{FQki} \\[1ex]
\mathrm{Qh} &= \mathrm{QC} \cdot \mathrm{FQh} \\[1ex]
\mathrm{Qha} &= \mathrm{Qh} - \mathrm{Qgu} \\[1ex]
\mathrm{Qlu} &= \mathrm{QC} \cdot \mathrm{FQlu} \\[1ex]
\mathrm{Qre} &= \mathrm{QC} \cdot \mathrm{FQre} \\[1ex]
\mathrm{Qpo} &= \mathrm{Qgu} \\[1ex]
\mathrm{GFR}_{\mathrm{base}} &= f_{\mathrm{GFR\_Qki}} \cdot \mathrm{Qki} \\[1ex]
\mathrm{Vlumen} &= 0.9 \cdot \mathrm{Vgu} \\[1ex]
\mathrm{Vki}_{\mathrm{plasma}} &= \mathrm{Vki} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vki}_{\mathrm{tissue}} &= \mathrm{Vki} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{plasma}} &= \mathrm{Vli} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vli}_{\mathrm{tissue}} &= \mathrm{Vli} \cdot \mathopen{}\left(1 - f_{\mathrm{tissue\_loss}}\right) \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{plasma}} &= \mathrm{Vlu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vlu}_{\mathrm{tissue}} &= \mathrm{Vlu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{plasma}} &= \mathrm{Vgu} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vgu}_{\mathrm{tissue}} &= \mathrm{Vgu} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{plasma}} &= \mathrm{Vre} \cdot \mathrm{Fblood} \cdot \mathopen{}\left(1 - \mathrm{HCT}\right) \\[1ex]
\mathrm{Vre}_{\mathrm{tissue}} &= \mathrm{Vre} \cdot \mathopen{}\left(1 - \mathrm{Fblood}\right) \\[1ex]
\mathrm{Ki}_{\mathrm{emp}} &= \frac{0.693}{\mathrm{ti}_{\mathrm{emp}}} \cdot 60 \\[1ex]
\mathrm{KI}_{\mathrm{\_GFR}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{GFR}_{\mathrm{base}} \\[1ex]
\mathrm{Aurine}_{\mathrm{emptot}} &= \mathrm{Aurine}_{\mathrm{emp}} + \mathrm{Aurine}_{\mathrm{eg}} \\[1ex]
\mathrm{KI}_{\mathrm{\_pg}} &= \mathrm{KI}_{\mathrm{\_fpg}} + \mathrm{KI}_{\mathrm{\_PG\_delta}} \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG\_fpg}} &= \mathrm{KI}_{\mathrm{\_RTG\_base}} + \mathrm{KI}_{\mathrm{\_RTG\_m\_fpg}} \cdot \mathopen{}\left(\mathrm{KI}_{\mathrm{\_fpg}} - \mathrm{KI}_{\mathrm{\_fpg\_healthy}}\right) \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG\_delta}} &= \mathrm{KI}_{\mathrm{\_RTG\_fpg}} \cdot \mathrm{KI}_{\mathrm{\_RTG\_max\_inhibition}} \\[1ex]
\mathrm{KI}_{\mathrm{\_RTG}} &= \mathrm{KI}_{\mathrm{\_RTG\_fpg}} - \frac{\mathrm{KI}_{\mathrm{\_RTG\_delta}} \cdot \mathrm{Cki}_{\mathrm{plasma\_emp}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}}}{\mathrm{KI}_{\mathrm{\_RTG\_E50}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}} + \mathrm{Cki}_{\mathrm{plasma\_emp}}^{\mathrm{KI}_{\mathrm{\_RTG\_gamma}}}} \\[1ex]
\mathrm{KI}_{\mathrm{\_UGE}} &= \frac{\mathrm{KI}_{\mathrm{\_glc\_urine}} \cdot \mathrm{KI}_{\mathrm{\_Mr\_glc}}}{\mathrm{KI}_{\mathrm{\_cf\_mg\_per\_g}}} \\[1ex]
\mathrm{GU}_{\mathrm{\_f\_absorption}} &= 1 + \mathrm{GU}_{\mathrm{\_fed}} \cdot \mathopen{}\left(\mathrm{GU}_{\mathrm{\_FED\_f\_absorption}} - 1\right) \\[1ex]
\mathrm{GU}_{\mathrm{\_absorption}} &= \mathrm{GU}_{\mathrm{\_f\_absorption}} \cdot \mathrm{GU}_{\mathrm{\_EMPABS\_k}} \cdot \mathrm{Vlumen} \cdot \mathrm{Cgu}_{\mathrm{emp}}
\end{aligned}
$$

## Reactions

| Rate | Id | Name | Equation | Modifiers |
| --- | --- | --- | --- | --- |
| $v_{\mathrm{transport\_lu\_emp}}$ | `transport_lu_emp` | transport empagliflozin | $\mathrm{Clu}_{\mathrm{plasma\_emp}} \rightleftharpoons \mathrm{Clu}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{transport\_re\_emp}}$ | `transport_re_emp` | transport empagliflozin | $\mathrm{Cre}_{\mathrm{plasma\_emp}} \rightleftharpoons \mathrm{Cre}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{transport\_re\_deep\_emp}}$ | `transport_re_deep_emp` | transport empagliflozin (deep tissue) | $\mathrm{Cre}_{\mathrm{plasma\_emp}} \rightleftharpoons \mathrm{Cre}_{\mathrm{deep\_emp}}$ |  |
| $v_{\mathrm{iv\_emp}}$ | `iv_emp` | iv empagliflozin | $\varnothing \longrightarrow \mathrm{Cve}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_emp}}$ | `Flow_ar_ki_emp` | inflow kidney empagliflozin | $\mathrm{Car}_{\mathrm{emp}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_emp}}$ | `Flow_ki_ve_emp` | outflow kidney empagliflozin | $\mathrm{Cki}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Cve}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_emp}}$ | `Flow_arli_li_emp` | arterial inflow liver empagliflozin | $\mathrm{Car}_{\mathrm{emp}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_emp}}$ | `Flow_arli_hv_emp` | flow arterial shunts | $\mathrm{Car}_{\mathrm{emp}} \longrightarrow \mathrm{Chv}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_emp}}$ | `Flow_po_li_emp` | outflow po empagliflozin | $\mathrm{Cpo}_{\mathrm{emp}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_emp}}$ | `Flow_po_hv_emp` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{emp}} \longrightarrow \mathrm{Chv}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_emp}}$ | `Flow_li_hv_emp` | outflow liver empagliflozin | $\mathrm{Cli}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Chv}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_emp}}$ | `Flow_hv_ve_emp` | outflow hepatic vein empagliflozin | $\mathrm{Chv}_{\mathrm{emp}} \longrightarrow \mathrm{Cve}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_emp}}$ | `Flow_ve_lu_emp` | inflow lung empagliflozin | $\mathrm{Cve}_{\mathrm{emp}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_emp}}$ | `Flow_lu_ar_emp` | outflow lung empagliflozin | $\mathrm{Clu}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Car}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_emp}}$ | `Flow_ar_gu_emp` | inflow gut empagliflozin | $\mathrm{Car}_{\mathrm{emp}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_emp}}$ | `Flow_gu_po_emp` | outflow gut empagliflozin | $\mathrm{Cgu}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Cpo}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_emp}}$ | `Flow_ar_re_emp` | inflow rest empagliflozin | $\mathrm{Car}_{\mathrm{emp}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_emp}}$ | `Flow_re_ve_emp` | outflow rest empagliflozin | $\mathrm{Cre}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Cve}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{Flow\_ar\_ki\_eg}}$ | `Flow_ar_ki_eg` | inflow kidney empagliflozin-glucuronide | $\mathrm{Car}_{\mathrm{eg}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_ki\_ve\_eg}}$ | `Flow_ki_ve_eg` | outflow kidney empagliflozin-glucuronide | $\mathrm{Cki}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Cve}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_arli\_li\_eg}}$ | `Flow_arli_li_eg` | arterial inflow liver empagliflozin-glucuronide | $\mathrm{Car}_{\mathrm{eg}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_arli\_hv\_eg}}$ | `Flow_arli_hv_eg` | flow arterial shunts | $\mathrm{Car}_{\mathrm{eg}} \longrightarrow \mathrm{Chv}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_po\_li\_eg}}$ | `Flow_po_li_eg` | outflow po empagliflozin-glucuronide | $\mathrm{Cpo}_{\mathrm{eg}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_po\_hv\_eg}}$ | `Flow_po_hv_eg` | flow portal shunts | $\mathrm{Cpo}_{\mathrm{eg}} \longrightarrow \mathrm{Chv}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_li\_hv\_eg}}$ | `Flow_li_hv_eg` | outflow liver empagliflozin-glucuronide | $\mathrm{Cli}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Chv}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_hv\_ve\_eg}}$ | `Flow_hv_ve_eg` | outflow hepatic vein empagliflozin-glucuronide | $\mathrm{Chv}_{\mathrm{eg}} \longrightarrow \mathrm{Cve}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_ve\_lu\_eg}}$ | `Flow_ve_lu_eg` | inflow lung empagliflozin-glucuronide | $\mathrm{Cve}_{\mathrm{eg}} \longrightarrow \mathrm{Clu}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_lu\_ar\_eg}}$ | `Flow_lu_ar_eg` | outflow lung empagliflozin-glucuronide | $\mathrm{Clu}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Car}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_ar\_gu\_eg}}$ | `Flow_ar_gu_eg` | inflow gut empagliflozin-glucuronide | $\mathrm{Car}_{\mathrm{eg}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_gu\_po\_eg}}$ | `Flow_gu_po_eg` | outflow gut empagliflozin-glucuronide | $\mathrm{Cgu}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Cpo}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{Flow\_ar\_re\_eg}}$ | `Flow_ar_re_eg` | inflow rest empagliflozin-glucuronide | $\mathrm{Car}_{\mathrm{eg}} \longrightarrow \mathrm{Cre}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{Flow\_re\_ve\_eg}}$ | `Flow_re_ve_eg` | outflow rest empagliflozin-glucuronide | $\mathrm{Cre}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Cve}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{KI\_\_EMP2EG}}$ | `KI__EMP2EG` | empagliflozin glucuronidation (EMP2EG) | $\mathrm{Cki}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Cki}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{KI\_\_EMPEX}}$ | `KI__EMPEX` | empagliflozin excretion (EMPEX) | $\mathrm{Cki}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Aurine}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{KI\_\_EGEX}}$ | `KI__EGEX` | EG excretion (EGEX) | $\mathrm{Cki}_{\mathrm{plasma\_eg}} \longrightarrow \mathrm{Aurine}_{\mathrm{eg}}$ |  |
| $v_{\mathrm{KI\_\_GLCEX}}$ | `KI__GLCEX` | glucose excretion (GLCEX) | $\mathrm{KI}_{\mathrm{\_fpg}} \longrightarrow \mathrm{KI}_{\mathrm{\_glc\_urine}}$ | $\mathrm{Cki}_{\mathrm{plasma\_emp}}$ |
| $v_{\mathrm{LI\_\_EMP2EG}}$ | `LI__EMP2EG` | empagliflozin glucuronidation (EMP2EG) | $\mathrm{Cli}_{\mathrm{plasma\_emp}} \longrightarrow \mathrm{Cli}_{\mathrm{plasma\_eg}}$ |  |
| $v_{\mathrm{GU\_\_EMPABS}}$ | `GU__EMPABS` | absorption empagliflozin (EMPABS) | $\mathrm{Cgu}_{\mathrm{emp}} \longrightarrow \mathrm{Cgu}_{\mathrm{plasma\_emp}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC}}$ | `GU__EMPEXC` | empagliflozin to the fecal transit (EMPEXC) | $\mathrm{Cgu}_{\mathrm{emp}} \longrightarrow \mathrm{GU}_{\mathrm{\_emp\_intestine\_0}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC\_0}}$ | `GU__EMPEXC_0` | fecal transit empagliflozin 0 | $\mathrm{GU}_{\mathrm{\_emp\_intestine\_0}} \longrightarrow \mathrm{GU}_{\mathrm{\_emp\_intestine\_1}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC\_1}}$ | `GU__EMPEXC_1` | fecal transit empagliflozin 1 | $\mathrm{GU}_{\mathrm{\_emp\_intestine\_1}} \longrightarrow \mathrm{GU}_{\mathrm{\_emp\_intestine\_2}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC\_2}}$ | `GU__EMPEXC_2` | fecal transit empagliflozin 2 | $\mathrm{GU}_{\mathrm{\_emp\_intestine\_2}} \longrightarrow \mathrm{GU}_{\mathrm{\_emp\_intestine\_3}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC\_3}}$ | `GU__EMPEXC_3` | fecal transit empagliflozin 3 | $\mathrm{GU}_{\mathrm{\_emp\_intestine\_3}} \longrightarrow \mathrm{GU}_{\mathrm{\_emp\_intestine\_4}}$ |  |
| $v_{\mathrm{GU\_\_EMPEXC\_4}}$ | `GU__EMPEXC_4` | fecal transit empagliflozin 4 | $\mathrm{GU}_{\mathrm{\_emp\_intestine\_4}} \longrightarrow \mathrm{Afeces}_{\mathrm{emp}}$ |  |
| $v_{\mathrm{GU\_\_dissolution\_emp}}$ | `GU__dissolution_emp` | dissolution empagliflozin | $\mathrm{GU}_{\mathrm{\_emp\_stomach}} \longrightarrow \mathrm{Cgu}_{\mathrm{emp}}$ |  |

The rates of the reactions are:

$$
\begin{aligned}
v_{\mathrm{transport\_lu\_emp}} &= \frac{\mathrm{ftissue}_{\mathrm{emp}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Clu}_{\mathrm{plasma\_emp}} \cdot \mathrm{Kp}_{\mathrm{emp}} - \mathrm{Clu}_{\mathrm{emp}}\right) \\[1ex]
v_{\mathrm{transport\_re\_emp}} &= \frac{\mathrm{ftissue}_{\mathrm{emp}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{plasma\_emp}} \cdot \mathrm{Kp}_{\mathrm{emp}} - \mathrm{Cre}_{\mathrm{emp}}\right) \\[1ex]
v_{\mathrm{transport\_re\_deep\_emp}} &= \frac{\mathrm{ftissue}_{\mathrm{deep\_emp}} \cdot \mathrm{BW}}{\mathrm{BW}_{\mathrm{ref}}} \cdot \mathopen{}\left(\mathrm{Cre}_{\mathrm{plasma\_emp}} \cdot \mathrm{Kp}_{\mathrm{deep\_emp}} - \mathrm{Cre}_{\mathrm{deep\_emp}}\right) \\[1ex]
v_{\mathrm{iv\_emp}} &= \frac{\mathrm{Ki}_{\mathrm{emp}} \cdot \mathrm{IVDOSE}_{\mathrm{emp}}}{\mathrm{Mr}_{\mathrm{emp}}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_emp}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_emp}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_emp}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_emp}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_emp}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_emp}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_emp}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_emp}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_emp}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_emp}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_emp}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_emp}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_emp}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{emp}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_emp}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{Flow\_ar\_ki\_eg}} &= \mathrm{Qki} \cdot \mathrm{Car}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_ki\_ve\_eg}} &= \mathrm{Qki} \cdot \mathrm{Cki}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{Flow\_arli\_li\_eg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_arli\_hv\_eg}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qha} \cdot \mathrm{Car}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_po\_li\_eg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_po\_hv\_eg}} &= f_{\mathrm{shunts}} \cdot \mathrm{Qpo} \cdot \mathrm{Cpo}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_li\_hv\_eg}} &= \mathopen{}\left(1 - f_{\mathrm{shunts}}\right) \cdot \mathopen{}\left(\mathrm{Qpo} + \mathrm{Qha}\right) \cdot \mathrm{Cli}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{Flow\_hv\_ve\_eg}} &= \mathrm{Qh} \cdot \mathrm{Chv}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_ve\_lu\_eg}} &= \mathrm{Qlu} \cdot \mathrm{Cve}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_lu\_ar\_eg}} &= \mathrm{Qlu} \cdot \mathrm{Clu}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{Flow\_ar\_gu\_eg}} &= \mathrm{Qgu} \cdot \mathrm{Car}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_gu\_po\_eg}} &= \mathrm{Qgu} \cdot \mathrm{Cgu}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{Flow\_ar\_re\_eg}} &= \mathrm{Qre} \cdot \mathrm{Car}_{\mathrm{eg}} \\[1ex]
v_{\mathrm{Flow\_re\_ve\_eg}} &= \mathrm{Qre} \cdot \mathrm{Cre}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{KI\_\_EMP2EG}} &= \mathrm{KI}_{\mathrm{\_f\_renal\_function}} \cdot \mathrm{KI}_{\mathrm{\_f\_oat3}} \cdot f_{\mathrm{ugt}} \cdot \mathrm{EMP2EG}_{\mathrm{k\_renal}} \cdot \mathrm{Vki}_{\mathrm{tissue}} \cdot \mathrm{Cki}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{KI\_\_EMPEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_fu\_emp}} \cdot \mathopen{}\left(1 + \mathrm{KI}_{\mathrm{\_f\_oat3}} \cdot \mathrm{KI}_{\mathrm{\_EMPSEC\_f}}\right) \cdot \mathrm{Cki}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{KI\_\_EGEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_EGEX\_f}} \cdot \mathrm{Cki}_{\mathrm{plasma\_eg}} \\[1ex]
v_{\mathrm{KI\_\_GLCEX}} &= \mathrm{KI}_{\mathrm{\_GFR}} \cdot 0.001 \cdot \mathrm{KI}_{\mathrm{\_RTG\_splay}} \cdot \ln\mathopen{}\left(1 + \exp\mathopen{}\left(\frac{\mathrm{KI}_{\mathrm{\_pg}} - \mathrm{KI}_{\mathrm{\_RTG}}}{\mathrm{KI}_{\mathrm{\_RTG\_splay}}}\right)\right) \\[1ex]
v_{\mathrm{LI\_\_EMP2EG}} &= \mathrm{LI}_{\mathrm{\_f\_oatp}} \cdot f_{\mathrm{ugt}} \cdot \mathrm{EMP2EG}_{\mathrm{k}} \cdot \mathrm{Vli}_{\mathrm{tissue}} \cdot \mathrm{Cli}_{\mathrm{plasma\_emp}} \\[1ex]
v_{\mathrm{GU\_\_EMPABS}} &= \mathrm{GU}_{\mathrm{\_F\_emp\_abs}} \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC}} &= \mathopen{}\left(1 - \mathrm{GU}_{\mathrm{\_F\_emp\_abs}}\right) \cdot \mathrm{GU}_{\mathrm{\_absorption}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC\_0}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vintestine\_0}} \cdot \mathrm{GU}_{\mathrm{\_emp\_intestine\_0}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC\_1}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vintestine\_1}} \cdot \mathrm{GU}_{\mathrm{\_emp\_intestine\_1}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC\_2}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vintestine\_2}} \cdot \mathrm{GU}_{\mathrm{\_emp\_intestine\_2}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC\_3}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vintestine\_3}} \cdot \mathrm{GU}_{\mathrm{\_emp\_intestine\_3}} \\[1ex]
v_{\mathrm{GU\_\_EMPEXC\_4}} &= \mathrm{GU}_{\mathrm{\_METEXC\_k}} \cdot \mathrm{GU}_{\mathrm{\_Vintestine\_4}} \cdot \mathrm{GU}_{\mathrm{\_emp\_intestine\_4}} \\[1ex]
v_{\mathrm{GU\_\_dissolution\_emp}} &= \frac{\frac{\mathrm{GU}_{\mathrm{\_Ka\_dis\_emp}}}{60} \cdot \mathrm{PODOSE}_{\mathrm{emp}}}{\mathrm{Mr}_{\mathrm{emp}}}
\end{aligned}
$$

## ODE system

The states change in time with the rates of the reactions and the rate rules:

$$
\begin{aligned}
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_emp}} - v_{\mathrm{Flow\_ki\_ve\_emp}} - v_{\mathrm{KI\_\_EMP2EG}} - v_{\mathrm{KI\_\_EMPEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_emp}} + v_{\mathrm{Flow\_po\_li\_emp}} - v_{\mathrm{Flow\_li\_hv\_emp}} - v_{\mathrm{LI\_\_EMP2EG}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_emp}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_lu\_emp}} + v_{\mathrm{Flow\_ve\_lu\_emp}} - v_{\mathrm{Flow\_lu\_ar\_emp}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_emp}} - v_{\mathrm{Flow\_gu\_po\_emp}} + v_{\mathrm{GU\_\_EMPABS}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_emp}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{transport\_re\_emp}} - v_{\mathrm{transport\_re\_deep\_emp}} + v_{\mathrm{Flow\_ar\_re\_emp}} - v_{\mathrm{Flow\_re\_ve\_emp}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_emp}} - v_{\mathrm{Flow\_arli\_li\_emp}} - v_{\mathrm{Flow\_arli\_hv\_emp}} + v_{\mathrm{Flow\_lu\_ar\_emp}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_emp}} - v_{\mathrm{Flow\_ar\_re\_emp}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Vve}} \Bigl( v_{\mathrm{iv\_emp}} + v_{\mathrm{Flow\_ki\_ve\_emp}} + v_{\mathrm{Flow\_hv\_ve\_emp}} - v_{\mathrm{Flow\_ve\_lu\_emp}} \\
&\quad {} + v_{\mathrm{Flow\_re\_ve\_emp}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_emp}} - v_{\mathrm{Flow\_po\_hv\_emp}} + v_{\mathrm{Flow\_gu\_po\_emp}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_emp}} + v_{\mathrm{Flow\_po\_hv\_emp}} + v_{\mathrm{Flow\_li\_hv\_emp}} - v_{\mathrm{Flow\_hv\_ve\_emp}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_lu\_emp}}}{\mathrm{Vlu}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_emp}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{deep\_emp}}}{\mathrm{d} t} &= \frac{v_{\mathrm{transport\_re\_deep\_emp}}}{\mathrm{Vre}_{\mathrm{tissue}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{emp}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_EMPEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cki}_{\mathrm{plasma\_eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_ki\_eg}} - v_{\mathrm{Flow\_ki\_ve\_eg}} + v_{\mathrm{KI\_\_EMP2EG}} - v_{\mathrm{KI\_\_EGEX}}}{\mathrm{Vki}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cli}_{\mathrm{plasma\_eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_li\_eg}} + v_{\mathrm{Flow\_po\_li\_eg}} - v_{\mathrm{Flow\_li\_hv\_eg}} + v_{\mathrm{LI\_\_EMP2EG}}}{\mathrm{Vli}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Clu}_{\mathrm{plasma\_eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ve\_lu\_eg}} - v_{\mathrm{Flow\_lu\_ar\_eg}}}{\mathrm{Vlu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{plasma\_eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_gu\_eg}} - v_{\mathrm{Flow\_gu\_po\_eg}}}{\mathrm{Vgu}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cre}_{\mathrm{plasma\_eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ar\_re\_eg}} - v_{\mathrm{Flow\_re\_ve\_eg}}}{\mathrm{Vre}_{\mathrm{plasma}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Car}_{\mathrm{eg}}}{\mathrm{d} t} &= \frac{1}{\mathrm{Var}} \Bigl( -v_{\mathrm{Flow\_ar\_ki\_eg}} - v_{\mathrm{Flow\_arli\_li\_eg}} - v_{\mathrm{Flow\_arli\_hv\_eg}} + v_{\mathrm{Flow\_lu\_ar\_eg}} \\
&\quad {} - v_{\mathrm{Flow\_ar\_gu\_eg}} - v_{\mathrm{Flow\_ar\_re\_eg}} \Bigr) \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cve}_{\mathrm{eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_ki\_ve\_eg}} + v_{\mathrm{Flow\_hv\_ve\_eg}} - v_{\mathrm{Flow\_ve\_lu\_eg}} + v_{\mathrm{Flow\_re\_ve\_eg}}}{\mathrm{Vve}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cpo}_{\mathrm{eg}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{Flow\_po\_li\_eg}} - v_{\mathrm{Flow\_po\_hv\_eg}} + v_{\mathrm{Flow\_gu\_po\_eg}}}{\mathrm{Vpo}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Chv}_{\mathrm{eg}}}{\mathrm{d} t} &= \frac{v_{\mathrm{Flow\_arli\_hv\_eg}} + v_{\mathrm{Flow\_po\_hv\_eg}} + v_{\mathrm{Flow\_li\_hv\_eg}} - v_{\mathrm{Flow\_hv\_ve\_eg}}}{\mathrm{Vhv}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Aurine}_{\mathrm{eg}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_EGEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Cgu}_{\mathrm{emp}}}{\mathrm{d} t} &= \frac{-v_{\mathrm{GU\_\_EMPABS}} - v_{\mathrm{GU\_\_EMPEXC}} + v_{\mathrm{GU\_\_dissolution\_emp}}}{\mathrm{Vlumen}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{Afeces}_{\mathrm{emp}}}{\mathrm{d} t} &= v_{\mathrm{GU\_\_EMPEXC\_4}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{KI}_{\mathrm{\_glc\_urine}}}{\mathrm{d} t} &= v_{\mathrm{KI\_\_GLCEX}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_emp\_intestine\_0}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_EMPEXC}} - v_{\mathrm{GU\_\_EMPEXC\_0}}}{\mathrm{GU}_{\mathrm{\_Vintestine\_0}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_emp\_intestine\_1}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_EMPEXC\_0}} - v_{\mathrm{GU\_\_EMPEXC\_1}}}{\mathrm{GU}_{\mathrm{\_Vintestine\_1}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_emp\_intestine\_2}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_EMPEXC\_1}} - v_{\mathrm{GU\_\_EMPEXC\_2}}}{\mathrm{GU}_{\mathrm{\_Vintestine\_2}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_emp\_intestine\_3}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_EMPEXC\_2}} - v_{\mathrm{GU\_\_EMPEXC\_3}}}{\mathrm{GU}_{\mathrm{\_Vintestine\_3}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{GU}_{\mathrm{\_emp\_intestine\_4}}}{\mathrm{d} t} &= \frac{v_{\mathrm{GU\_\_EMPEXC\_3}} - v_{\mathrm{GU\_\_EMPEXC\_4}}}{\mathrm{GU}_{\mathrm{\_Vintestine\_4}}} \\[1ex]
\frac{\mathrm{d} \,\mathrm{PODOSE}_{\mathrm{emp}}}{\mathrm{d} t} &= -v_{\mathrm{GU\_\_dissolution\_emp}} \cdot \mathrm{Mr}_{\mathrm{emp}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{IVDOSE}_{\mathrm{emp}}}{\mathrm{d} t} &= -v_{\mathrm{iv\_emp}} \cdot \mathrm{Mr}_{\mathrm{emp}} + \mathrm{Ri}_{\mathrm{emp}} \qquad \text{(rate rule)} \\[1ex]
\frac{\mathrm{d} \,\mathrm{cum}_{\mathrm{dose\_emp}}}{\mathrm{d} t} &= \mathrm{Ri}_{\mathrm{emp}} \qquad \text{(rate rule)}
\end{aligned}
$$
