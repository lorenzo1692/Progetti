# PFC — Poloidal Field Coils

Pipeline di design delle Poloidal Field Coils (PFC): port modulare di
`PF_opt_VNS.m` sul modello di TF e CS (vedi il manuale PFC). "VNS" è il nome
del caso studio, non un algoritmo: la ricerca è un'enumerazione esaustiva
con filtri di fattibilità.

## Come si usa

```matlab
addpath(genpath('MADE-2.0/matlab/PFC'))   % oppure: cd PFC; main_WP_PFC_design
main_WP_PFC_design
```

Il main chiede, in ordine:

1. il file dei parametri (`input/WP_PFC_input_template.xlsx`, da copiare e
   modificare: tutti i parametri sono nella colonna *Value*);
2. il file di geometria e Ampere-turns di tutta la macchina (righe CS, PF1..PF6,
   plasma; colonne R, Z, dr, dz e una per scenario, es.
   `Baseline_VNS_07_2026_V3_CREATE.xlsx`);
3. il file delle forze assiali per scenario (`FZ_PF.xlsx`), oppure Invio per
   calcolarle dalle correnti di scenario (`physics/compute_scenario_forces.m`).

Per ogni bobina PF1..PF6 e ogni `WP_h` dello sweep produce la tabella dei design
point fattibili (`<tag>_PF<n>_results_<data>.xlsx`), poi permette di sfogliarli,
sceglierne uno, salvarlo e verificarlo con il FEM.

## Struttura

- `input/` — template Excel con i 57 parametri (geometria, punto di design,
  numerica, materiali, **fatica FCGR**, **FEM**, check geometrici).
- `io/` — lettura dell'input (`read_machine_input`), geometria/scenario
  (`read_coil_geometry`), forze assiali (`read_axial_force`).
- `physics/` — dimensionamento conduttore (`size_conductor_cicc`) e jacket
  (`size_cicc_cable`), campo e forze (`emag_field_forces`,
  `estimate_peak_field`), sforzo ad anello (`eqv_stress_coil_ring_cicc`),
  fatica (`fcgr`, `ParametriY`), matrice di induttanza
  (`compute_coupling_matrix`), forze di scenario (`compute_scenario_forces`),
  materiali in `physics/materials/`.
- `search/` — `generate_combinations` e `scan_wp_designs`.
- `postprocess/` — `browse_solutions`, `plot_solution`.
- `fem/` — verifica FE assialsimmetrica della soluzione scelta
  (`fem_pfc_verify`), vedi sotto.

Dipendenze condivise nella cartella superiore (`MADE-2.0/matlab/`): `xbr.m`,
`xbz.m`, `xlm.m`, `heat_balance_cicc_ode.m` (il main aggiunge il path).

## Fatica FCGR come opzione di dimensionamento

`fcgr_mode` nel file di input: `0` = spenta, `1` = calcola `plasma_cycles`
(migliaia di cicli) per ogni design accettato, `2` = il jacket viene fatto
crescere finché `plasma_cycles >= plasma_cycles_min`. Le costanti del materiale
(`fcgr_C0`, `fcgr_m`, ...) sono nel file di input (316LN e JK2LB indicati nelle
descrizioni).

## Verifica FEM (`fem/`)

`fem_pfc_verify(riga_DATA, p, geom)` ricostruisce il winding pack del design point
e lo risolve con un FE assialsimmetrico (r,z), quadrilateri bilineari, materiale
ortotropo omogeneizzato, per:

- il caso di design (bobina sola, corrente `Iop`, come ipotizza il modello
  analitico);
- ogni scenario di plasma: corrente della bobina scalata sugli Ampere-turns di
  scenario e campo di tutte le altre bobine (CS, PF, plasma).

I carichi sono le forze di Lorentz `J x B` nei punti di Gauss. Il confronto stampa
sforzo di hoop, verticale, Tresca, forze radiali/assiali contro il modello
analitico (`eqv_stress_coil_ring_cicc`, `emag_field_forces`).
`fem_selftest` verifica il solver contro la soluzione di Lamé (errore < 0.1%),
l'equilibrio e una colonna soggetta al peso proprio.

Stato: v1. Materiale omogeneizzato (niente modello del singolo turn), nessuna
struttura di supporto, vincolo verticale a scelta (`fem_bc`).
