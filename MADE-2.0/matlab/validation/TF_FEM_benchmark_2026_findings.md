# TF FEM benchmark 2026 — validation findings

Source: ANSYS APDL 2D generalized-plane-strain run, 21-Sep-2026 (`TFBM_*.txt`
export), geometry patched per `patches/change_notes.json` from
`TF_FEM_benchmark_v0_1`. 100 turns, 2 field grades, 11 physical rows (6 rows
of 10 turns at 32.6 mm cell height + 5 rows of 8 turns at 22.8 mm), 12 TF
coils, Iop = 62.3 kA, JT = 3.5 mm, cable fillet r_SC = 4 mm, cable modulus
E_cbl = 10 GPa, TREF=293 K / TUNIF=4.2 K, frictional WP-case contacts.

## Confirmed FEM peak stress (nodal-averaged, global Cartesian axes)

| Component | Node | SINT | SEQV | Discretization bound (SINT) |
|---|---|---:|---:|---:|
| Jacket (mat 2) | 183613 | **979.6 MPa** | 898.0 MPa | up to ~1095 MPa |
| Case/Vault (mat 10) | 217646 | **773.6 MPa** | 698.4 MPa | up to ~792 MPa |

These match the historical plot values (980 / 774 MPa) almost exactly,
confirming the earlier geometry-index bugs (`H_GRADES`, `Rj_` closure) were
metadata/summary bugs that did not change the actual solved mesh.

The historical SCL section (nodes 183612/183605) sits inside **row 6, the
last row of grade 1** (Y ≈ 0.958–0.962 m, inside the 0.9569–0.9895 m range of
row 6), i.e. at the **grade 1 / grade 2 transition**, not in the innermost
grade as the current analytical model assumes (see below).

## Forward evaluation of the current analytical formulas

`validate_TF_FEM_benchmark_2026.m` evaluates `size_cicc_cable.m`'s stiffness
formula and `size_case_vault.m`'s stress-balance formula at this exact FEM
geometry (no sizing search), for two material assumptions.

**Round 1 — checking only the last layer (original code):**

| Test | E_cbl | r_SC | S_T_JT (Jacket) | S_T_VT (Vault) |
|---|---:|---:|---:|---:|
| A — FEM-matched | 10 GPa | 4 mm | 450.2 MPa | 671.8 MPa |
| B — current tool defaults | 0.1 GPa | 5 mm | 487.4 MPa | 695.7 MPa |
| FEM reference (SINT) | — | — | 979.6 MPa | 773.6 MPa |

**Round 2 — checking every layer, still no transition correction:** 458.5 MPa
(grade 1 rows) vs 450.2 MPa (grade 2 rows) for Test A — checking every layer
instead of only the last one closes only **~2%** of the gap. The FEM peak
location (row 6, the grade1/grade2 boundary) is confirmed but a uniform
per-layer formula cannot reproduce it: geometry is identical within each
grade, so the formula gives essentially the same answer everywhere in that
grade.

**Round 3 — with `SCF_transition_provisional` = 3.15 applied to the two
layers adjacent to a grade change (rows 6 and 7 here):**

| Test | S_T_JT (Jacket) | vs FEM SINT (979.6 MPa) |
|---|---:|---:|
| A — FEM-matched (calibration case) | **979.0 MPa** | −0.06% |
| B — current tool defaults | **1068.8 MPa** | +9.1% |

## Reading the gap

1. **Vault/Case stress was already in the right ballpark** (−10 to −13% vs
   FEM SINT, unaffected by this change) — the stiffness-ratio coupling in
   `size_case_vault.m` is not badly wrong.
2. **The Jacket gap is not a location bug.** Checking every layer instead of
   only the last barely moves the number (2%). The FEM peak comes from local
   bending/deformation at the grade-transition stiffness discontinuity,
   which a lumped stiffness-network model cannot produce by construction.
3. **The E_cbl mismatch (0.1 vs 10 GPa) is a minor factor on its own**
   (~8% swing) but **compounds with the transition SCF**: calibrated under
   FEM-matched materials (Test A, exact match), the same SCF overshot by
   +9% when combined with the tool's old E_cbl_LTS=0.1 GPa default
   (Test B). This was a second, independent piece of evidence (after the
   review's own audit) that `E_cbl_LTS` should be corrected to match the
   real cable modulus - **now done**: the input template default is 10 GPa
   (see "Implemented" below), so Test A is now simply "current tool
   defaults" and the +9% overshoot no longer applies.
4. `SCF_transition_provisional` (now in `input/WP_TF_input_template.xlsx`,
   category "Structural corrections (provisional)") is a **single-point
   calibration** - it makes the tool stop grossly under-predicting the
   Jacket peak, but it is not validated across geometries, currents, or
   grade counts. Treat it as a stopgap, not physics.

## Update — all hardcoded parameters exposed, r_SC now a clamp rule

Every remaining hardcoded numeric literal in the physics/search modules
(sizing-loop steps, geometric feasibility thresholds, the cable corner
fillet radius, the vacuum-vessel margin factor) is now an input-template
parameter — see the "Numerical settings" and "WP dimensioning" categories.
The cable corner fillet radius `r_SC` is no longer a flat constant: it now
follows `r_SC = clamp(JT, r_SC_min, r_SC_max)` (default 2-6 mm), i.e. it
tracks the jacket thickness up to 6 mm and is fixed beyond that. Re-running
the forward evaluation with this rule (JT=3.5mm -> r_SC=3.5mm here) and the
corrected E_cbl_LTS=10GPa gives **S_T_JT=978.5 MPa, S_T_VT=669.7 MPa**,
essentially unchanged from the round-3 numbers above (r_SC has a small
effect, as already suggested by the E_cbl-driven Test A/B comparison).

## Implemented (this pass)

- `search/scan_wp_designs.m`: the "Primary radial stress (Pm+Pb)" check now
  loops over **every** layer (previously only `var = n_layers`) and keeps
  the worst case; layers adjacent to a grade transition (from `jump_grade`)
  are multiplied by `p.SCF_transition_provisional`.
- New parameter `SCF_transition_provisional` (default 3.15) added to the
  input template and to `read_machine_input.m`'s required list.
- `input/WP_TF_input_template.xlsx`: `E_cbl_LTS` default corrected from
  0.1 GPa to 10 GPa, matching the FEM's actual cable structural modulus
  (ANSYS material 1, EX=1e10 Pa) confirmed in `TFBM_model.txt`. Test A/B
  above are no longer "FEM-matched vs tool defaults" - both point to the
  same value now.

## Suggested next step

Phase 2 (agreed, not yet started): replace `SCF_transition_provisional`
with a physics-based local-bending correction at the layer-stiffness
discontinuity (e.g. treating adjacent grades as elastically-coupled rings/
beams with a compatibility condition at the interface), so the correction
generalizes to geometries other than this one benchmark case. Re-run this
same forward evaluation after any change, and revisit the E_cbl_LTS default
question above.

Not yet addressed by this comparison (still open, lower priority for now):
thermal contraction (TUNIF=4.2K cooldown), contact nonlinearity, and the
WP_h/Rj_ summary formula in `search/scan_wp_designs.m` not including the
10×0.5mm inter-row insulation gaps that the layer-by-layer Re/Ri chain does
account for (~5 mm effect on Rj_/DTF, i.e. a few percent of the 135.4 mm
nose thickness in this case).
