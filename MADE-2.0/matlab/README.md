# MADE 2.0 — dimensionamento del winding pack TF (MATLAB)

Punto di ingresso: **`main_WP_TF_design.m`**. Tutti i parametri della
macchina sono nel file di input Excel (`input/WP_TF_input_template.xlsx`,
da copiare per ogni macchina): il codice non contiene parametri di
macchina. Il manuale completo (funzioni, teoria, validazione) è
`docs/manuale_web.html`.

## Struttura

| Cartella | Contenuto | File principali |
|---|---|---|
| `input/` | File di input Excel | `WP_TF_input_template.xlsx` |
| `io/` | Lettura dell'input | `read_machine_input.m` |
| `physics/` | Grandezze di macchina, geometria dei turn, campo 2D, dimensionamento, modello FE meccanico | `compute_operating_params.m`, `wp_envelope.m`, `wp_turn_geometry.m` (Rect / RIS), `wp_field_at_points.m`, `compute_discrete_field_profile.m`, `wp_peak_field_fast.m`, `wp_field_calibration.m` + `wp_field_calibrated.m` (campo di picco tarato), `size_cicc_cable.m`, `size_case_vault.m`, `wp_mech_surrogate.m` |
| `conductor/` | Conduttore CICC: scelta del superconduttore, strand, rame, hot spot | `cicc.m`, `cicc_params.m`, `heat_balance_cicc_ode.m`, `Ic_Nb3Sn.m`, `Ic_NbTi.m`, `Ic_sst33.m`, `jc_ybco.m` |
| `search/` | Scansione dello spazio di progetto; ridimensionamento di una soluzione; dimensionamento accoppiato con il modello globale 3D | `generate_combinations.m`, `scan_wp_designs.m`, `resize_design_point.m`, `couple_global_sizing.m` |
| `postprocess/` | Scelta della soluzione, plot, ricarica, export | `browse_solutions.m`, `plot_*.m`, `load_design_point.m`, `export_ansys_input.m`, `ask_save_plots.m`, `plot_wp_section_manual_example.m` (configurazione a mano) |
| `coil3d/` | Modulo 3D: shape bending-free, campo 3D Biot–Savart, export; modello globale a travi e shell di tutto il sistema TF (controparte di STR_360) costruito dalla soluzione scelta | `tf3d_from_design.m`, `tf_bending_free_shape.m`, `biot_savart_segments.m`, `tf3d_global_from_design.m`, `tf3d_global_model.m`, `tf3d_global_solve.m` |
| `validation/` | Confronti con ANSYS e verifiche riproducibili; `results/` riferimenti e log, `tools/` script Python di estrazione | `validate_mech_surrogate_2026.m`, `verify_ris_rect_fem.m`, `validate_tf3d.m` |
| `docs/` | Manuale, report, note di consegna | `manuale_web.html`, `RIS_E_RESIDUO_FEM.md`, `HANDOFF*.md` |
| `legacy/` | Codice della versione precedente, **non usato** e fuori dal path | vedi `legacy/README.md` |

## Uso rapido

```matlab
main_WP_TF_design                                   % scansione completa, scelta, plot, FE, 3D, modello globale e dimensionamento accoppiato, export
[row, p] = load_design_point('<tag>_results_<data>.xlsx', idx);   % ricarica una soluzione salvata
out = wp_mech_surrogate(row, p); plot_wp_mech_surrogate(out, p);  % verifica meccanica FE
R = verify_ris_rect_fem(p, 'verify_out');           % verifica RIS / Rect
```
