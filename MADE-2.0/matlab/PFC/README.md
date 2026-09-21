# PFC — Poloidal Field Coils

Scaffolding per la pipeline di design delle Poloidal Field Coils (PFC),
speculare alla pipeline del TF (vedi `main_WP_TF_design.m` e le cartelle
`input/`, `io/`, `physics/`, `search/`, `postprocess/` in
`MADE-2.0/matlab/`):

- `input/` — file di input (es. template Excel) con i parametri delle PFC.
- `io/` — lettura e validazione dell'input in una struct di parametri.
- `physics/` — calcolo del punto operativo e dimensionamento del WP/CICC.
- `search/` — generazione e scan delle combinazioni di design candidate.
- `postprocess/` — analisi, tabelle e grafici delle soluzioni trovate.

Da popolare seguendo lo stesso pattern modulare del dominio TF.
