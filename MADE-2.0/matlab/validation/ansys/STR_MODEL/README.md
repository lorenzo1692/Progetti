# 2D model of the TF inner leg (ANSYS APDL)

Electromagnetic and structural 2D model of one wedged TF coil (inner leg),
built from the design parameters exported by MADE
(`export_ansys_input.m` -> `input/ParametriTF.f`).

## How to run

Launch `RUN.dat` in batch from this folder, for example:

    ansys2xx -b -i RUN.dat -o RUN.out

`RUN.dat` creates the subfolder `output\` (`/MKDIR`) and makes it the
working directory (`/CWD`, all MPI processes). Every file the run writes goes there:

| File | Content |
|---|---|
| `EM_2D.db`, `EM_2D.rst` | EM model and solution (Lorentz forces) |
| `STR_2D.db`, `STR_2D.rst` | structural model and solution |
| `*.png` | plots of `Solution_EM.f` and `postpro/Plot.f` |
| `CHECK_TF.txt` | resultant of the Lorentz forces and T_bf/2 [MN] |
| `TFBM_*.txt` | diagnostic exports (`EXPORT_TF_*.mac`) |

Solver scratch files are deleted at the end. The batch log (`RUN.out`)
stays where ANSYS was launched.

The model files are read through the parameter `MDIR` (`..\` seen from
`output\`). A model file started by hand from this folder sets
`MDIR = '.\'` and writes to the current directory, as before.

`RUN_BENCHMARK.dat` re-opens an existing solution in `output\` and
writes only the `TFBM_*.txt` exports.

## Flow

| File | Role |
|---|---|
| `RUN.dat` | batch driver: EM model + solution, structural model + solution, cleanup |
| `EM_MODEL.f` | EM geometry: conductors and air |
| `Solution_EM.f` | PLANE233 solution, current per turn, EM export and plots |
| `STR_MODEL.f` | structural geometry: conductors, jackets, insulation, fillers, ground and wedge insulation, case |
| `BC/Contacts.f` | element types and contact pairs |
| `Solution_STR.f` | flank constraints, vertical force (generalized plane strain), Lorentz forces (LDREAD), cool-down, solution |
| `postpro/Plot.f` | plots, insulation fatigue index |
| `input/ParametriTF.f` | design parameters (from MADE) |
| `input/Materiali.f` | materials at 4.2 K |
| `geom/*_trpz.f` | WP with layers of different width (the current MADE designs) |
| `geom/*_RIS.f` | round-in-square WP (`WP_TOPOLOGY = 'RIS_only'`) |

The other `geom/` macros (equal-width layers, pancakes, buffer) are not
called by the current flow; each one says so in its header.

## Changes (October 2026)

1. **Jacket thickness varying between layers.** The cable half-width in
   `geom/TF_cab_trpz.f`, `TF_jck_trpz.f`, `TF_tins_trpz.f` and
   `TF_fill_trpz.f` used the jacket of layer 1 (`JT_w(1,1)`), while the
   jacket wall used the jacket of its own layer. Layers with a thicker
   jacket than layer 1 overflowed their cell and overlapped the next turn,
   and ANSYS failed. All four macros now use `JT_w(n,k)`. Also fixed:
   - `R_N = R(n,k)` replaces `R(j,k)`;
   - the jacket internal radius is now set per grade (`CFRG`, written by
     `export_ansys_input`, as in MADE: clamp(JT, r_SC_min, r_SC_max)).
   Parameter files without `CFRG` keep the single `CFR`.
2. **Output folder** `output\`, as described above. The hard-coded
   `/CWD` to the network folder in `STR_MODEL.f` and `EM_MODEL.f` was
   removed: the model now runs wherever it is launched.
3. **Comments** in English, with a header in every file (purpose,
   caller, inputs and outputs). Typos and misleading section titles were
   fixed, for example "DOUBLE PANCAKE INSULATION" in the inter-layer
   insulation macros. No command was changed apart from points 1, 2 and 4.
4. **Air mesh of the EM model** (`geom/TF_air.f`). The last node merge,
   after the air is mirrored, used the default tolerance of 0.1 mm. In a
   design with variable jacket it merged two distinct nodes of a small air
   triangle and collapsed it ("aspect ratio 1.E+20", "zero or negative
   determinant of the Jacobian"): 6 errors and the batch run stopped. It
   now merges with 1e-6 m, as the final merges in `TF_fill_trpz.f` and
   `TF_cfill.f`: the nodes that must merge (symmetry axis, mirrored air on
   the cables of the other half) coincide exactly.

## Open points

- `STR_MODEL.f` calls `geom/TF_gins_panck.f` for `WP_TOPOLOGY = 'RIS_only'`,
  but the file is not in `geom/`: restore it before running a RIS model.
- `geom/nuovo 1` is an unnamed scratch copy of an inter-layer insulation
  macro and is not used.
- Output folder: created with `/MKDIR` and entered with `/CWD`. The first
  version used `/SYS,mkdir`, which fails when the model sits on a network
  (UNC) path, because cmd.exe cannot start there. If `/CWD` fails anyway,
  `RUN.dat` detects it (`/INQUIRE` + `STRPOS`) and runs in the model
  folder as before.
