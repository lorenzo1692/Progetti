# Extraction of the FEM reference values (validation of wp_mech_surrogate)

The reference numbers written in `../validate_mech_surrogate_2026.m` come from
the ANSYS exports `TFBM_*.txt` (nodal-averaged, material-restricted stresses,
global Cartesian axes) through these scripts (Python 3, numpy, scipy,
matplotlib). They are kept here so the extraction can be repeated on a new
FEM run.

| Script | What it does |
|---|---|
| `parse_fem.py <dir> <tag>` | Reads nodes and elements from `TFBM_global.txt` -> `<tag>_mesh.pkl` |
| `fem_lin.py` | Parsers of the stress/displacement listings, Tresca, linearization along a line of FEM nodes |
| `cmp_sections.py <tag> <dir> <sections.csv>` | Linearizes the FEM jacket stresses on exactly the sections exported by the surrogate (`out.sections`, with the layer appended as last column), by interpolation on the jacket elements, and compares Pm and Pm+Pb |
| `turn_peaks.py <tag> <dir> <surrogate.mat> "<geometry dict>" <JT>` | Maximum SINT of the FEM jacket per turn and per layer, against the surrogate peaks |
| `case_scl.py` | Linearization of the FEM case stresses along given lines (the surrogate case SCLs, `out.case_scl(q).P0/P1`) |
| `fem_lin_generic.py` | Linearization on the straight jacket walls at 1/4, 1/2, 3/4 of each wall from the design geometry |

Typical sequence for a new FEM run:

```
python3 parse_fem.py EM_2D007 em7
# run wp_mech_surrogate on the same geometry in MATLAB, then export
#   dlmwrite('sections.csv', [out.sections, [out.turn(out.sections(:,1)).layer]'], 'precision', '%.8e');
python3 cmp_sections.py em7 EM_2D007 sections.csv
python3 turn_peaks.py em7 EM_2D007 surr.mat "dict(Ri_=...,dps=0.02,GIT=0.005,INS=0.0005,n_turns=[...],Cond_w=[...],Cond_h=[...])" 0.003
```
