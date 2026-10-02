# CS — Central Solenoid

Pipeline di design del Central Solenoid: port modulare di `CS_opt_VNS.m` sul
modello di TF e PFC (vedi il Manuale CS). "VNS" è il nome del caso studio, non un
algoritmo: la ricerca è un'enumerazione esaustiva con filtri di fattibilità.

## Come si usa

```matlab
cd MADE-2.0/matlab/CS
main_WP_CS_design
```

Il main chiede il file dei parametri (`input/WP_CS_input_template.xlsx`, da
copiare e modificare) e il file scenario della macchina (es.
`Baseline_VNS_22_07_2025.xlsx`: colonne R, Z, dr, dz e una per scenario). Poi scansiona
turns/layers/Iop, salva i design fattibili, permette di sceglierne uno, produce i plot
(scatter di scan, sezione di un modulo e del CS completo colorata per |B|) e, a
richiesta, lo verifica con il FEM dello stack.

## Struttura

- `input/` — template Excel (52 parametri, incluse le opzioni FEM).
- `io/` — `read_machine_input`, `read_scenario_currents` (bound Ampere-turns),
  `read_coil_geometry` (per il FEM).
- `physics/` — `size_conductor_cicc`, `size_cicc_cable`, `eqv_stress_coil_cicc`,
  `emag_field_forces`, `fcgr`, `ParametriY`, `heat_balance_cicc_ode`, materiali in
  `physics/materials/`.
- `search/` — `generate_combinations`, `scan_wp_designs`.
- `postprocess/` — `browse_solutions`, `plot_solution`, `plot_cs_section`.
- `fem/` — verifica assialsimmetrica dello stack (vedi sotto).
- `test/` — `cs_smoke_test` (scan senza dialoghi, anche in Octave) e shim Octave.

Dipendenze condivise nella cartella superiore (`MADE-2.0/matlab/`): `xbr.m`, `xbz.m`,
`xlm.m` (il main aggiunge il path).

## Forza assiale

`Fz_MN` (colonna dei risultati) è la compressione calcolata per ogni candidato da
`emag_field_forces` con tutto lo stack (`full_stack=1`): somma delle Fz della metà
inferiore, cioè la forza che attraversa il piano medio. Non è più il valore fisso 35 MN
del legacy.

## Verifica FEM dello stack (`fem/`)

`fem_cs_verify(riga_DATA, p, g, geom, opts)` costruisce lo stack di `n_moduli` moduli
(isolamento di massa, regione dei turns omogeneizzata ortotropa, piastre distanziatrici
in acciaio fra i moduli), calcola le forze di Lorentz `J x B` nei punti di Gauss e
risolve il FE assialsimmetrico (r,z). Casi di carico: progetto (stack solo, tutti i
moduli a `Iop`) e ogni scenario (correnti di modulo, PF e plasma da `geom`). Il confronto
stampa, per modulo, forze, hoop, verticale e Tresca contro `eqv_stress_coil_cicc` e
`emag_field_forces`, e la compressione attraverso ogni piastra. Supporto: inerzia (`fem_bc=1`)
o base fissa con precarico in testa (`fem_bc=2`, `fem_preload_MN`).
`fem_cs_selftest(p, g)` controlla equilibri e coerenza con `emag_field_forces`.

Modulo dettagliato: con `fem_detail_module = m` il modulo `m` è modellato turn per turn
(isolamento di turn, parete del jacket, cavo; cella rettangolare, raggio d'angolo del
jacket trascurato) mentre gli altri restano omogenei. Le sue tensioni di hoop,
verticale e Tresca sono lette direttamente dagli elementi del jacket. Sul design di
prova (modulo 3, caso di progetto) il hoop del jacket è circa il 2-3% sotto il recupero
omogeneo, il Tresca circa il 9% sotto: la versione omogenea è prudente.

Esportazione per il modello ANSYS (`CS_model_2`): `export_apdl_params(riga, p, g, dir, tag)`
scrive `<tag>_Parametri_CS_design.lgw` (override del blocco CS di `Parametri_TCM_*.lgw`) e
`<tag>_DESIGN.csv` (scenario di progetto, formato di `STR/input/SN.csv`: tutti i moduli a
`Iop`, niente PF né plasma). Non ancora provato in ANSYS.
