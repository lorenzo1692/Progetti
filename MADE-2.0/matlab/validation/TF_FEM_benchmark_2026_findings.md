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
formula and `size_case_vault.m`'s stress-balance formula **once**, at this
exact FEM geometry (no sizing search), for two material assumptions:

| Test | E_cbl | r_SC | S_T_JT (Jacket) | S_T_VT (Vault) |
|---|---:|---:|---:|---:|
| A — FEM-matched | 10 GPa | 4 mm | **450.2 MPa** | **671.8 MPa** |
| B — current tool defaults | 0.1 GPa | 5 mm | **487.4 MPa** | **695.7 MPa** |
| FEM reference (SINT) | — | — | 979.6 MPa | 773.6 MPa |

## Reading the gap

1. **Vault/Case stress is in the right ballpark** (−10 to −13% vs FEM SINT).
   The stiffness-ratio coupling in `size_case_vault.m` is not badly wrong.
2. **Jacket stress is under-predicted by roughly a factor of 2** (−50 to
   −54%). This is the dominant error, not the Case side.
3. **The E_cbl mismatch (0.1 vs 10 GPa) is a minor factor** (8% difference
   between Test A and B), contrary to the initial hypothesis raised in
   chat — the dominant gap is elsewhere.
4. The most likely explanation: `size_cicc_cable.m`'s "Primary radial
   stress (Pm+Pb)" check only evaluates the **last** layer (innermost grade),
   while the FEM's true peak sits at the **grade transition** (row 6/7
   boundary) — a location this formula never checks. A lumped
   membrane+bending-via-stiffness-ratio model also cannot capture a true 2D
   bending/contact stress concentration at a material discontinuity like
   that transition.

## Suggested next step

Before building a fully coupled Jacket/Case equilibrium solver, first
address the Jacket-side gap, since it is the larger of the two:
- Evaluate the "Primary radial stress" check at **every** layer (not only
  the last), and in particular at grade transitions, keeping the maximum.
- Consider whether the lumped Pm+Pb formula needs a stress-concentration or
  local-bending correction near grade transitions, informed by this FEM
  case.
- Re-run this same forward evaluation after any change, to track whether
  the Jacket gap actually closes.

Not yet addressed by this comparison (still open, lower priority for now):
thermal contraction (TUNIF=4.2K cooldown), contact nonlinearity, and the
WP_h/Rj_ summary formula in `search/scan_wp_designs.m` not including the
10×0.5mm inter-row insulation gaps that the layer-by-layer Re/Ri chain does
account for (~5 mm effect on Rj_/DTF, i.e. a few percent of the 135.4 mm
nose thickness in this case).
