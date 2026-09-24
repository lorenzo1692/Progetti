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

## Discrete field model vs FEM (peak BSUM on conductor)

`physics/compute_discrete_field_profile.m` on the benchmark geometry
(11 layers, n_turns 6x10 + 5x8, Cond_w 46 mm, Cond_h 32.6/22.8 mm, JT 3.5 mm,
Iop 62.3 kA, n_TF 12):

| Model | Peak B on conductor |
|---|---|
| FEM (BSUM) | 13.489 T |
| Smeared Ampere (scan, B_TF) | 13.489 T |
| Old: line filaments at centroids, no self-field | 12.92 T (-4.2%) |
| New: coil-1 turns as uniform rectangles (self-field incl.), peak over 3x3 points per cable | 13.480 T (-0.07%) |

The old model also under-predicted the low-field layers by up to ~50%
(layer 11: 1.7 T vs 3.4 T), because there the turn's own self-field
(~mu0*I/(2*pi*a) ~ 1 T) dominates.

## Design 7 (EM_2D007 FEM run, 24-Sep-2026)

Scan design point 7: 13 layers x 8 turns (104 turns), 3 grades (4 + 3 + 6
layers, Cond_h 34.2 / 26.1 / 24.2 mm, Cond_w 43.0 mm, JT 3.1/3.0/3.0 mm),
Iop = 64.628 kA, Ri_ = 1.2594 m, Rk_ = 0.7301 m. Script:
`validate_TF_FEM_design7_2026.m` (uses `forward_eval_wp_stress.m`).

**Geometry notes.** (1) The FEM nose is 133.0 mm, the scan row says 139.0 mm:
the scan's `Rj_` omitted the (n_layers-1) inter-layer insulation gaps
(12 x 0.5 mm) - **fixed** in `search/scan_wp_designs.m`. The FEM (built with
`export_ansys_input`'s WPH) already had the right stack. (2) The FEM mesh has
JT = 3.0 mm in grade 1 (cable 35.0 x 26.2 mm) where the row has 3.1 mm.

### Field - discrete model confirmed, smeared model under-predicts

| | Peak B on conductor |
|---|---|
| FEM BSUM | **14.482 T** |
| Discrete model (`compute_discrete_field_profile`) | **14.472 T** (-0.07%) |
| Smeared B_TF used by the scan | 13.483 T (-6.9%) |

Field each grade's cable was **sized at** vs the real peak on that grade:
13.48 -> 14.47 T (+7%), 9.33 -> 11.24 T (+20%), 6.22 -> 8.97 T (+44%).
With `Ic_Nb3Sn` at the real field, the three grades carry only
**0.64 / 0.60 / 0.53 x Iop** (sized for 1.00 x Iop): the smeared
(Ampere, continuous-shell) field is not safe for sizing the cable,
especially in the low-field grades. (The first benchmark matched only
because its B_TF was taken from the FEM.) Note: `TFBM_em.txt` has no B
data ("requested B data is not available") - the EM results were not
loaded before PRNSOL; the peak was read from EM_2D002.png.

### Stress

| | Jacket | Case |
|---|---:|---:|
| FEM SINT | **1124.4 MPa** (L13) | **800.5 MPa** |
| Scan row | 651.1 (-42%) | 666.2 (-17%) |
| Analytical, FEM geometry, B smeared, transition SCF | 647 | 669 |
| Analytical, FEM geometry, B discrete, transition SCF | 717 (-36%) | 742 (**-7%**) |
| Analytical, FEM geometry, B discrete, no SCF | 355 | 742 |

- **Case:** with the correct field, -7% vs FEM (benchmark: -13%). The vault
  model is acceptable; most of the scan's -17% came from the field.
- **Jacket: the transition SCF is not the right model.** FEM SINT per
  layer is high everywhere and rises toward the nose within each grade:
  674, 864, 932, **977** (last of grade 1), 879, 929, **958**, 939, 965,
  984, 1004, 1025, **1124** MPa (L13, outermost turn, cable fillet corner
  next to the nose/side wall). The first benchmark shows the same pattern
  (613 -> 980 in grade 1, 695 -> 876 in grade 2): there the peak happened
  to be the last layer of grade 1, which is why a "transition" factor fitted
  it. The FEM/analytical ratio (no SCF) goes from 1.9 (L1) to 3.2 (L13),
  so 3.15 at the transitions gives the right number only by coincidence.
- **What the formula misses** (FEM jacket components, central turns):
  axial SZ 160-290 MPa (analytical S_z = 187 MPa - OK); radial compression
  in the side walls that **accumulates with depth**, +45 -> -312 MPa (the
  analytical radial term is the same in every layer); **toroidal wedging
  compression** -200 to -305 MPa in the top/bottom walls (absent from the
  analytical jacket formula); corner shear/bending at the fillet (SXY up to
  500 MPa at the peak node).
- A one-parameter depth-dependent fit (accumulated radial pressure) over
  both FEM runs still leaves ~11% rms / 28% max error: a proper per-layer
  jacket model (accumulated radial + wedging + corner factor) needs to be
  developed and calibrated on both runs (24 layer data points).

## 2D FE mechanical surrogate (physics/wp_mech_surrogate.m)

The lumped analytical model cannot give the jacket stress: the FEM jacket
stress comes from the true equilibrium of the section (radial load
accumulating towards the nose, toroidal wedging compression, bending of
the walls around the fillets, WP/case and cable/jacket sliding). A light
FE model that solves that equilibrium from the design point alone was
built and validated against both ANSYS runs:

- geometry as the FEM: rounded turns (cable fillet r_SC, jacket outer
  radius r_SC+JT, insulation r_SC+JT+tins, corner filler), inter-layer and
  ground insulation, case with flat plasma side, radial flanks and nose arc
  of radius Rk_/cos(pi/n_TF), cavity = convex hull of the layers + GIT,
  wedge insulation;
- boundary conditions as the FEM: zero normal displacement on the flanks,
  free sliding (the FEM wedge insulation is constrained only normally, so
  its mu = 0.2 carries no load), generalized plane strain with Fz = T_bf;
- contacts as the FEM: WP/case and cable/jacket unilateral with Coulomb
  friction mu = 0.2 (penalty + stick/slip iterations);
- loads: Lorentz force from the discrete field model (per-turn force
  within 0.03% of the ANSYS LDREAD loads), cool-down 293 -> 4.2 K with the
  benchmark CTEs (orthotropic insulation), axial force.

Validation (`validate_mech_surrogate_2026.m`, surrogate/FEM):

| | Benchmark (11 layers) | Design 7 (13 layers) |
|---|---|---|
| axial strain eps_z | +0.4% | +0.0% |
| jacket straight-wall sections, Pm / Pm+Pb (mean +- std) | 1.00+-0.04 / 1.01+-0.04 (1400) | 0.99+-0.06 / 0.99+-0.07 (1456) |
| jacket fillet sections (45 deg), Pm / Pm+Pb | 0.98 / 0.98 | 0.95 / 0.94 |
| jacket peak in the fillet, per layer | 0.92 .. 1.05 | 0.89 .. 1.05 |
| global jacket peak | 923 vs 980 MPa (-6%) | 1020 vs 1124 MPa (-9%) |
| case nose / vault SCL Pm | -3% / -2% | -5% / -3% |
| case side-wall SCLs | -2 .. +2% | +3% |

Peaks are located in the same fillets as in the FEM. What matters most,
learned while building it: the rounded corners (turns touch only along the
flat parts: with square corners the WP is 20-30% too stiff toroidally) and
the cable/jacket frictional contact (bonded: peaks -5..-26%; frictionless:
edge-turn peaks +27%).

**Figure of merit (after the code review, docs/RISPOSTA_REVIEW.md).**
Linearization alone does not classify stresses, so the surrogate now
solves two load cases: primary P (Lorentz + axial force, no cool-down) and
total P+Q (all loads, the cool-down being the secondary load). On the
linearized Tresca stresses of every jacket wall and fillet section and of
the case SCLs it checks Pm(P) <= Sm, (Pm+Pb)(P) <= 1.5 Sm and
(Pm+Pb)(P+Q) <= 3 Sm, with Sm_jacket / Sm_case explicit in the Excel input
(default 667 MPa, TO BE CONFIRMED against the design code). The fillet
peak (P+Q) is reported for information (local stress: fatigue / FEM
check). Every result carries validity checks (mesh area, Jacobian, contact
convergence per load step, residual, global and axial equilibrium, SCL
coverage); a failed check marks the figure of merit INVALID.

Design 7: jacket Pm(P) 618 MPa (0.93 Sm), (Pm+Pb)(P) 1107 MPa
(1.11 x 1.5 Sm, NOT satisfied: without the cool-down clamping the jacket
walls carry the Lorentz load in bending), (Pm+Pb)(P+Q) 882 MPa
(0.44 x 3 Sm), peak 1020 MPa; case Pm 623 MPa (0.93 Sm). With contacts the
P / P+Q split is approximate (the primary case is solved with its own
contact state), which is conservative. The earlier "criteria satisfied"
verdict did not classify the cool-down and is superseded.

Benchmark: jacket Pm(P) 595 MPa (0.89 Sm), (Pm+Pb)(P) 981 MPa
(0.98 x 1.5 Sm), (Pm+Pb)(P+Q) 796 MPa (0.40 x 3 Sm); case Pm 660 MPa
(0.99 Sm) - satisfied.

The validation script now has explicit acceptance bands (eps_z +-2%;
per-layer peak [0.85 1.10], Pm [0.90 1.10], Pm+Pb [0.85 1.10]; case SCLs
[0.90 1.05] nose/vault, [0.95 1.05] side walls, [0.95 1.20] plasma-side
plate), requires the validity checks to pass and ends with an error on any
failure. Run time: ~4 min per design in Octave with both load cases,
expected well under a minute in MATLAB; used on the chosen design point
(main_WP_TF_design step 5b), not inside the combinatorial scan.
