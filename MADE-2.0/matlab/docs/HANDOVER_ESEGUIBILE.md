# Handover: eseguibile di MADE 2.0 (MATLAB)

Documento per la chat che deve **costruire un eseguibile standalone**
basato sulla versione attuale di MADE 2.0. Leggerlo tutto prima di
modificare il codice.

## 0. Obiettivo e vincoli

**Obiettivo:** un'applicazione eseguibile, compilata con MATLAB Compiler
(`.exe` su Windows, binario su Linux, con MATLAB Runtime), che faccia lo
stesso lavoro di `main_WP_TF_design.m`:
1. input da file Excel;
2. scansione del winding pack;
3. scelta della soluzione;
4. plot;
5. FE meccanico facoltativo;
6. modulo 3D facoltativo;
7. export ANSYS;
8. salvataggio delle figure.

**Vincoli, da rispettare alla lettera:**
- **La fisica non si tocca.** Tutto ciò che sta in `physics/`, `conductor/`, `search/` e `coil3d/` resta com'è. Stesso vale per soglie, default e tolleranze. L'eseguibile è un involucro, non una revisione del modello.
- **Stessi risultati dello script.** A parità di file di input, l'eseguibile deve dare la stessa tabella di soluzioni (§6). Ogni differenza è un bug dell'involucro.
- **`main_WP_TF_design.m` deve continuare a funzionare** come oggi, in modalità interattiva da MATLAB.
- **Non inventare risultati.** Se non hai MATLAB/Compiler per eseguire o compilare, dillo, e consegna codice e istruzioni di build da verificare.

## 1. Versione di partenza

- Repository `lorenzo1692/progetti`, branch `claude/code-improvement-3cmng0`, commit `b68732e` o successivi, che aggiungono solo documentazione (oppure lo ZIP `MADE-2.0-matlab.zip` di questa versione).
- Cartella del tool: `MADE-2.0/matlab/`. Mappa delle cartelle in `README.md`; manuale completo (teoria, funzioni, validazione) in `docs/manuale_web.html`.

```
matlab/
  main_WP_TF_design.m     punto di ingresso interattivo (da conservare)
  README.md               mappa delle cartelle
  input/                  WP_TF_input_template.xlsx (foglio "Input": Category, Name, Value, Unit, Description)
  io/                     read_machine_input.m
  physics/                campo 2D, taratura del campo, dimensionamento, FE meccanico
  conductor/              cicc.m, heat_balance_cicc_ode.m, leggi Ic
  search/                 generate_combinations.m, scan_wp_designs.m
  postprocess/            scelta soluzione, plot, ricarica, export ANSYS, salvataggio plot
  coil3d/                 modulo 3D (ricerca)
  validation/             confronti con ANSYS e verifiche (NON va nell'eseguibile)
  docs/                   manuale e note
  legacy/                 codice non usato, fuori dal path (NON va nell'eseguibile)
```

## 2. Il flusso attuale (`main_WP_TF_design.m`)

| Passo | Cosa fa | Funzioni | Interazione oggi | Output |
|---|---|---|---|---|
| 1 | Sceglie e fa confermare il file di input | — | `input()` × 2 | — |
| 2–4 | Legge l'input, parametri derivati, inviluppo, combinazioni, **scansione** (con taratura iniziale del campo se `field_model` = 1) | `read_machine_input`, `compute_operating_params`, `wp_envelope`, `generate_combinations`, `scan_wp_designs` | nessuna | `DATA` (table), `field_cal`; figura della taratura |
| 5 | Salva i risultati | `writetable`, `save` | nessuna | `<tag>_results_<data>.xlsx` + `.mat` (DATA, p, field_cal) in `pwd` |
| 5 | Scelta della soluzione | `browse_solutions` | `uifigure`/`uitable` + `ginput` + `input()` | `sel_idx` |
| 5 | Plot della soluzione | `plot_solution`, `plot_wp_section`, `plot_wp_section_bfield`, `plot_wp_diagnostics`, `plot_hotspot_transient` | nessuna | figure; `plot_solution` salva PNG/Excel della riga in `pwd` |
| 5b | FE meccanico (facoltativo) | `wp_mech_surrogate`, `plot_wp_mech_surrogate`, `plot_conductor_zoom` | `input()` | figure; se non valido, banner "NOT VALIDATED" |
| 6 | Modulo 3D (facoltativo) | `tf3d_from_design`, `export_tf3d_results`, `export_tf3d_shape`, `plot_tf3d` | `input()` | cartella `<tag>_tf3d_<data>` in `pwd`, figure |
| 7 | Export ANSYS (facoltativo) | `export_ansys_input` | `input()` | file APDL in `pwd` |
| 8 | Salvataggio figure | `ask_save_plots` | `input()` × 2 | PNG/FIG nella cartella scelta |

Tempi indicativi, misurati in Octave; MATLAB è più veloce:

| Operazione | Tempo in Octave |
|---|---|
| Scansione | ~4–8 s per candidato arrivato al dimensionamento (una scansione completa del template sono ore) |
| Taratura del campo | 5–6 min |
| FE meccanico | 3–6 min |
| Modulo 3D | ~35 s |

L'interfaccia deve quindi:
- mostrare l'avanzamento, oggi con una riga di console auto-sovrascritta in `scan_wp_designs`;
- non bloccarsi senza segnali.

## 3. Cosa impedisce oggi la compilazione e come risolverlo

1. **`input()` e `ginput`.** In un eseguibile senza console (e in ogni app GUI) non funzionano. Si trovano in:
   - `main_WP_TF_design.m`: righe con `input(` per file, conferma, FE, 3D, export;
   - `postprocess/browse_solutions.m`: `ginput` + `input`;
   - `postprocess/ask_save_plots.m`: `input` × 2.

   **Soluzione:** separare la logica dall'interazione (§4). Le funzioni interattive restano per lo script; l'eseguibile usa le versioni non interattive.
2. **`addpath` / `rmpath`.** Nel codice compilato sono vietate e danno errore: proteggerle con `if ~isdeployed`. Nell'eseguibile i file entrano con `mcc`/`compiler.build` (§5), e `legacy/` e `validation/` non vanno inclusi.
3. **Percorsi.**
   - `mfilename('fullpath')` e `this_dir` puntano dentro l'archivio compilato: per il template usare `fullfile(ctfroot, ...)` oppure includerlo con `-a` e trovarlo con `which`/`exist`.
   - Tutti gli output oggi vanno in `pwd` (`plot_solution`, `export_ansys_input`, `export_tf3d_*`, `ask_save_plots`, risultati del main). Nell'eseguibile serve una **cartella di output esplicita**, scelta dall'utente, passata a queste funzioni: quasi tutte hanno già l'argomento `out_dir`/`folder`.
4. **Figure.** In modalità batch crearle con `'Visible','off'` e salvarle su file. In GUI possono restare visibili. `savefig` (.fig) funziona anche compilato, ma PNG è sufficiente.
5. **`uifigure`/`uitable`** (browse_solutions) funzionano compilati; `ginput` no in una GUI: la selezione va fatta cliccando la riga della tabella o il punto dello scatter (callback).
6. **File di dati.** `conductor/Ic_sst33.m` carica file `.mat` solo per le opzioni angolari 1 e 3, che `cicc` non usa (usa `[3,4]`, analitiche). Nessun file di dati è quindi necessario. Il compilatore potrebbe però segnalare le dipendenze: ignorarle o gestirle, senza cambiare le opzioni usate.
7. **Toolbox.** Il codice usa solo MATLAB base:
   - `ode45`, `fminsearch`, `integral`, `besseli`;
   - `delaunayTriangulation`, `containers.Map`, `readtable`/`writetable`, `uifigure`/`uitable`.

   Verificarlo con `matlab.codetools.requiredFilesAndProducts('main_WP_TF_design.m')` prima di compilare; serve solo MATLAB Compiler.
8. **Warning di `ode45`** ("Solving was not successful…") dentro `cicc`: sono frequenti e noti. Non trasformarli in errori, e nella GUI non mostrarli uno per uno: vanno nel log.

## 4. Architettura proposta

Obiettivo: un nucleo non interattivo e due involucri.

```
made_run(cfg)                 NUCLEO, nessuna interazione, nessun addpath
   ├─ legge cfg.input_file, esegue passi 2–5 (scansione, salvataggi)
   ├─ selezione: cfg.select = indice | 'best' (criterio esplicito) | [] (solo scansione)
   ├─ passi facoltativi secondo cfg.run_fe, cfg.run_3d, cfg.export_ansys, cfg.save_plots
   ├─ scrive tutto in cfg.out_dir (+ log di testo)
   └─ restituisce res: DATA, field_cal, sel_idx, mech, tf3d, elenco dei file scritti
made_batch(input_file, out_dir, ...)   eseguibile da riga di comando (name/value)
made_app                               GUI (App Designer o uifigure programmatica)
main_WP_TF_design.m                    resta interattivo; può chiamare made_run
```

**GUI** (`made_app`), requisiti minimi:
1. **Input:** scelta del file Excel, con il template come default. Mostrare la tabella dei parametri (Category, Name, Value, Unit, Description), modificabile e salvabile in un nuovo file. Nessun parametro va scritto nel codice.
2. **Opzioni:** `field_model` (0/1/2) e `field_verify` sono già nel file Excel. Aggiungere i flag dei passi facoltativi: FE, 3D, export ANSYS, salvataggio figure.
3. **Esecuzione:** barra di avanzamento della scansione, con analizzati, accettati ed ETA, e un pulsante Annulla.
   - Per l'avanzamento e l'annullamento serve un callback facoltativo in `scan_wp_designs`: un solo argomento opzionale, senza cambiare il resto.
   - Un'intestazione di funzione con un argomento in più è l'unica modifica ammessa ai file della fisica.
4. **Risultati:**
   - tabella delle soluzioni, ordinabile, e scatter Iop–Rk con selezione per click;
   - poi le figure della soluzione in schede.
5. **FE:** mostrare esito dei controlli e figura di merito. Se `out.valid` è falso, avviso ben visibile come nei plot attuali ("NOT VALIDATED"): la figura di merito non va usata.
6. **Output:** cartella scelta dall'utente, contenente:
   - risultati `.xlsx`/`.mat`;
   - figure PNG;
   - export ANSYS e 3D;
   - un file di log con la versione (commit) e i parametri usati.

## 5. Build

Esempio, da adattare alla release MATLAB disponibile:

```matlab
% dalla cartella matlab/, con il path impostato (senza legacy/ e validation/)
opts = compiler.build.StandaloneApplicationOptions('made_app.m', ...
    'ExecutableName', 'MADE2', ...
    'AdditionalFiles', {'input/WP_TF_input_template.xlsx'}, ...
    'TreatInputsAsNumeric', 'off');
compiler.build.standaloneApplication(opts);
% versione batch
compiler.build.standaloneApplication('made_batch.m', 'ExecutableName', 'MADE2_batch', ...
    'AdditionalFiles', {'input/WP_TF_input_template.xlsx'});
```

Documentare:
- la release MATLAB e la versione di MATLAB Runtime necessarie;
- come si installa il Runtime;
- come si lancia la versione batch, con un esempio di riga di comando.

## 6. Verifiche di accettazione

Da fare in MATLAB, **prima** con lo script e **poi** con l'eseguibile, e da riportare con i numeri ottenuti:
1. **Validazione FE contro ANSYS** (solo script, non compilata): `validation/validate_mech_surrogate_2026.m` deve dare **PASSED**.
2. **Verifica RIS/Rect:** `R = verify_ris_rect_fem(p, 'verify_out')`, tutti i controlli superati.
   - Design 7 (Rect): Pm+Pb 882 MPa, picco 1020 MPa (tutti i carichi); Pm+Pb primario 1107 MPa → FoM non soddisfatta, valida.
   - RIS: Pm+Pb 571 MPa, picco 630 MPa, FoM soddisfatta.
3. **Scansione ridotta** (design 7): template con `field_model` = 0, 13 layer × 8 turn, spessore laterale `env.lateral_w_min + 6*p.lateral_w_step`.
   - Deve dare 1 soluzione con Iop 64.6 kA e grade a 13.48 / 9.334 / 6.223 T.
   - Il picco discreto sulla soluzione (`wp_peak_field_fast`) deve valere 14.482 T; ANSYS dà 14.482 T.
4. **Stessa scansione ridotta con `field_model` = 1:** confrontare `DATA` dello script e dell'eseguibile colonna per colonna. Devono coincidere, a meno dell'arrotondamento dell'ultima cifra.
5. **Esecuzione completa della GUI** sul template:
   - scansione, scelta, FE, 3D, export, salvataggio;
   - controllare che tutti i file finiscano nella cartella di output e nessuno in `pwd`.

## 7. Da NON fare in questo lavoro

Sono decisioni aperte, fuori dallo scopo dell'eseguibile:
- non correggere la formula di Ampère in `compute_operating_params` (`2*pi*RTFi - dr_plasma_side`), non legare `p_rs` al campo tarato, non toccare lo SCF;
- non cambiare `cicc`, i default del file di input o le soglie dei controlli del FE;
- non includere `legacy/` né `validation/` nell'eseguibile;
- non "sistemare" il modulo 3D: è marcato di ricerca e va lasciato com'è.

## 8. Riferimenti

- `README.md`: mappa delle cartelle.
- `docs/manuale_web.html`: manuale completo, compresi il flusso, il campo tarato (§6.9) e i controlli del FE (§8.9).
- `docs/RIS_E_RESIDUO_FEM.md`: conduttore RIS, residuo del FE, campo di picco, correzione di `cicc`.
- `legacy/README.md`: cosa c'è in `legacy/` e perché non si usa.
