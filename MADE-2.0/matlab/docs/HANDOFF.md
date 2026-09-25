# MADE 2.0 – documento di passaggio consegne

Questo file serve a chi riprende il lavoro sul codice (persona o assistente AI) senza aver seguito lo sviluppo. Contiene:

- cosa fa il codice;
- cosa è stato fatto e validato;
- le convenzioni da rispettare;
- i numeri di riferimento;
- le decisioni aperte;
- il prossimo lavoro: l'integrazione del passaggio al 3D, con l'analisi dei file originali già fatta.

Da leggere insieme a:

- `docs/REPORT_codice_e_modello_meccanico.md`: guida a tutte le funzioni e descrizione dettagliata del modello meccanico FE;
- `docs/RISPOSTA_REVIEW.md`: revisione del codice, punti C00–C08;
- `validation/TF_FEM_benchmark_2026_findings.md`: diario dei confronti con ANSYS, con i numeri.

---

## 1. Contesto

MADE 2.0 è un tool MATLAB per il **dimensionamento del winding pack (WP) delle bobine di campo toroidale (TF)** di un tokamak (progetto VNS 2026, 12 bobine, conduttori CICC Nb3Sn/NbTi/REBCO).

- Repository: `lorenzo1692/Progetti`, cartella `MADE-2.0/matlab/`.
- Branch di sviluppo: `claude/code-improvement-3cmng0`. Nessuna pull request aperta.
- Lingua: commenti nel codice in inglese, documentazione e dialogo con l'utente in italiano.
- L'utente è un ingegnere di magneti per la fusione: usa ANSYS per le verifiche FEM e MATLAB per il dimensionamento.

### Catena di calcolo

```
input/WP_TF_input_template.xlsx      (80 parametri, colonne Category/Name/Value/Unit)
  -> io/read_machine_input           -> struct p
  -> physics/compute_operating_params-> struct g (RTFi, RTFo, NI, k_bf, r_bf, B_PHI_TF, ...)
  -> physics/wp_envelope             -> env
  -> search/generate_combinations    -> combT
  -> search/scan_wp_designs          -> tabella DATA (una riga per soluzione fattibile)
  -> salvataggio .xlsx + .mat, browse_solutions -> sel_idx
  -> plot (sezione, campo, diagnostica, hot spot)
  -> passo 5b: physics/wp_mech_surrogate (FE 2D) + plot_wp_mech_surrogate
  -> passo 6:  postprocess/export_ansys_input (APDL)
```

Il punto d'ingresso è `main_WP_TF_design.m`. Una "soluzione" è una riga `row` di `DATA`, o una struct con gli stessi campi. È l'input di tutte le funzioni a valle.

**Campi principali di `row`:**

- scalari: `Iop`, `n_layers`, `Ri_`, `Rk_`, `Tau_discharge`, `S_T_VT`, `S_T_JT`;
- vettori per layer: `n_turns`, `Cond_w`, `Cond_h`, `JT`, `type_cable`, `N_Sc`, `N_Cu`, `THS`, `B_grade`.

### Convenzioni geometriche (identiche al FEM ANSYS)

- Sezione 2D della gamba interna: **x toroidale, y radiale**, asse macchina nell'origine, bobina 1 centrata sull'asse +y.
- Il layer 1 è quello lato plasma (raggio maggiore).
- `row.Ri_` è il raggio della faccia lato plasma del case della gamba interna.
- Cava del WP: parte da `Ri_ - dr_plasma_side`, poi isolante di massa `GoundIns` (il nome del parametro ha questo refuso, non va cambiato).
- Layer:
  ```
  Re(1) = Ri_ - dr_plasma_side - GoundIns
  Ri(k) = Re(k) - Cond_h(k)
  Re(k+1) = Ri(k) - INS_grades
  ```
- `row.Rk_` è il raggio interno del case (nose); l'arco del nose ha raggio `Rk_/cos(pi/n_TF)`.
- Un turn è una cella `Cond_w × Cond_h`. Dall'esterno verso l'interno: isolante di spira `turn_insulation_nominal`, jacket `JT`, cavo con raccordo `r_SC = clamp(JT, r_SC_min, r_SC_max)`.

---

## 2. Stato attuale (tutto committato e pushato)

L'ultimo commit con codice è `11315df` ("Address code review C00-C08 ...").

### 2.1 Campo magnetico discreto (validato)

- `physics/wp_field_at_points.m`:
  - campo 2D dei cavi della bobina analizzata come sezioni arrotondate a densità di corrente uniforme, scomposte esattamente in rettangoli con Biot–Savart analitico; l'auto-campo è incluso;
  - le altre bobine sono filamenti.
- `physics/compute_discrete_field_profile.m`: picco per turn su 48 punti del contorno più una griglia 3×3. Restituisce `x, y, layer, B_smooth, B_center, B_discrete, ripple, r_SC`.
- Risultati: **14.476 T contro 14.482 T del FEM** (design 7) e 13.488 contro 13.489 T (benchmark).

### 2.2 Modello meccanico FE 2D (validato)

`physics/wp_mech_surrogate.m` (circa 1400 righe, funzioni locali) è un FEM autonomo:

- **Formulazione:** deformazione piana generalizzata con ε_z globale e Fz = T_bf; settore 2π/n_TF.
- **Mesh:** Q8 curvi per jacket, isolanti e filler; T6 Delaunay per cavo, case e isolante di massa.
- **Contatti:** unilaterali con attrito di Coulomb incrementale, per le interfacce WP/case e cavo/jacket; i fianchi hanno solo il vincolo normale.
- **Carichi:** raffreddamento, Lorentz dal campo discreto, forza assiale.
- **Due casi di carico:** P (Lorentz + assiale) e P+Q (tutti i carichi).
- **Linearizzazione (SCL)** su ogni sezione del jacket e su 7 linee del case.
- **Criteri:** Pm ≤ Sm, Pm+Pb ≤ 1.5 Sm (sul caso P), P+Q ≤ 3 Sm.
- **Controlli di validità:** area, Jacobiano, convergenza dei contatti, residuo, equilibrio globale e assiale, copertura SCL. Se uno fallisce, la figura di merito è marcata INVALID.

**Validazione** (`validation/validate_mech_surrogate_2026.m`, bande di accettazione esplicite, esito PASS/FAIL): **PASSATA** su entrambi i casi ANSYS.

| Rapporto modello FE / FEM | Benchmark (11 layer) | Design 7 (13 layer) |
|---|---|---|
| εz | +0.3% | 0.0% |
| Picco jacket per layer | 0.92–1.05 | 0.89–1.05 |
| Pm per layer | 0.96–1.02 | 0.93–0.99 |
| Pm+Pb per layer | 0.90–1.04 | 0.855–1.02 (margine minimo: layer 1) |
| Nose / vault del case | 0.97 / 0.98 | 0.95 / 0.97 |

**Figura di merito del design 7** (Sm = 667 MPa, **da confermare** con l'utente e il codice di progetto):

| Criterio | Jacket | Esito |
|---|---|---|
| Pm primario | 618 MPa (0.93 Sm) | soddisfatto |
| Pm+Pb primario | **1107 MPa = 1.11 × 1.5 Sm** | **non soddisfatto** |
| P+Q | 882 MPa (0.44 × 3 Sm) | soddisfatto |

Il case è soddisfatto: Pm 623 MPa, 0.93 Sm.

Tempi: circa 4 minuti per design in Octave, con i due casi di carico.

### 2.3 Scansione

- Barra di avanzamento: soluzioni totali, analizzate, fattibili, rimanenti, tempo stimato.
- Correzioni già fatte:
  - raggio interno del WP (`Rj_`);
  - il vault ora è confrontato con `S_amm_VT` (punto C07).
- La scansione usa ancora **modelli analitici**: il campo spalmato (lineare per layer) e il fattore `SCF_transition_provisional = 3.15`, calibrato su un solo FEM.

### 2.4 Altro

- `postprocess/plot_hotspot_transient.m`: T(t) di hot spot per ogni grade.
- `postprocess/export_ansys_input.m`: file parametri APDL della soluzione.
- `validation/tools/*.py`: estrazione dei riferimenti dagli export ANSYS (vedi il README della cartella).

---

## 3. Come lavorare sul codice

- **Test:** eseguiti in **GNU Octave 8.4** (non c'è MATLAB nell'ambiente di sviluppo). Accorgimenti:
  - le funzioni locali negli *script* non sono visibili in Octave: usare file-funzione;
  - `readtable` non funziona: si usa una struct `p` generata dall'Excel con uno script Python (`openpyxl`);
  - `clim` → `set(gca,'CLim',...)`; `c.Label.String` → `ylabel(c,...)`; niente `sgtitle`;
  - `run()` cambia la cartella corrente: usare `addpath` con percorsi assoluti.
- **Il codice deve restare compatibile con MATLAB**, che è l'ambiente dell'utente.
- **Nuovi parametri:** vanno nell'Excel (`input/WP_TF_input_template.xlsx`), nella categoria giusta, con default letto tramite `get_or`/`isfield`, così i file di input vecchi continuano a funzionare.
- **Dopo ogni modifica a campo o meccanica:** rieseguire `validation/validate_mech_surrogate_2026.m` (circa 8 minuti in Octave). Deve finire con `PASSED`. Non modificare i riferimenti FEM né le bande per far passare un test.
- **Consegna all'utente:** a ogni gruppo di modifiche l'utente riceve lo zip della cartella `MADE-2.0/matlab`.
- **Commit:** su `claude/code-improvement-3cmng0`, messaggi descrittivi.

### Bug noto non ancora corretto

`physics/compute_operating_params.m` contiene:

```matlab
B_PHI_TF = g.Mu_0*p.n_TF*NI_MA*1e6/(2*pi*g.RTFi - p.dr_plasma_side);
```

Probabilmente mancano le parentesi: dovrebbe essere `2*pi*(g.RTFi - p.dr_plasma_side)`. L'effetto è piccolo (0.02 m contro 2π·1.2 m) ma cambia il campo spalmato della scansione. **Da confermare con l'utente prima di correggere**, perché sposta i risultati salvati.

---

## 4. Decisioni aperte con l'utente

1. **C06 – integrazione nella scansione.** Scegliere una o più strade:
   - (a) dimensionare ogni grade sul picco di campo discreto, iterando tra dimensionamento e campo;
   - (b) usare il modello FE come filtro obbligatorio sulle soluzioni fattibili;
   - (c) tarare sul modello FE una stima analitica veloce da usare al posto di `SCF_transition_provisional`.
2. **Sm:** confermare Sm = 667 MPa per jacket e case, e la classificazione del raffreddamento come carico secondario. Il verdetto sul design 7 dipende da entrambe.
3. **Bug `B_PHI_TF`** (vedi §3).

---

## 5. PROSSIMO LAVORO: passaggio al 3D a valle di MADE

### 5.1 Richiesta dell'utente

> "questi file vanno sistemati ma contengono il primo passaggio al 3D dopo aver scelto la sezione. mettili a valle di made"

I file originali sono in `coil3d/legacy/`, copiati senza modifiche. `ColorQuiver.m` è identico a quello già presente nella radice.

| File | Cosa fa |
|---|---|
| `TF_shape_opt_BF_Ncoils.m` | Script: forma D "bending free" per un numero finito di bobine, poi approssimazione con 3 archi tangenti, curve di offset di case e WP, export per CAD in mm |
| `bendingfree_opt.m` | Iterazione della forma bending free; a ogni iterazione chiama **ANSYS** (`MAIN.dat`, Biot–Savart con SOURC36) per il campo toroidale sulla linea della bobina |
| `MAIN.dat` | Macro APDL usata da `bendingfree_opt` |
| `emag_calculation_TFC.m` | Campo 3D con una bobina = un filamento sulla D analitica (Princeton D); Biot–Savart regolarizzato (Hurwitz); matrice di induttanza via potenziale vettore |
| `emag_calculation_TFC_n_wires.m` | Come sopra, con un pacchetto NL×NT di filamenti valutati sulla bobina 1; parametri DEMO |
| `emag_calculation_TFC_n_wires_full.m` | Pacchetto NL×NT su tutte le bobine e tentativo di induttanza completa |

### 5.2 Problemi trovati nei file originali (verificati leggendo il codice)

1. **Dipendenza da ANSYS.** `bendingfree_opt` lancia `C:\Program Files\ANSYS Inc\v221\...\ANSYS221.exe` a ogni iterazione e legge `B_.csv`. Funziona solo su quella macchina Windows.
2. **Forma costruita sulla linea sbagliata.** La forma bending free è calcolata per `r1 = Ri`, `r2 = Re`, cioè per la faccia lato plasma del case. Poi WP e case vengono ricavati per offset. La condizione di bending free vale invece per la **linea baricentrica della corrente** (centroide del WP).
3. **Sezione fittizia in ANSYS.** In `MAIN.dat` la sorgente ha sezione `WP_w = 1.0` × `WP_h = 0.5` m, fissa e non collegata al WP reale.
4. **Fit a 3 archi fragile** (`TF_shape_opt_BF_Ncoils.m`, blocco `while max_distance >= 0.40`):
   - le equazioni 7 e 8 hanno un errore di segno: `... - (xA - x)^2 + (yA - y)^2` invece di `- (yA - y)^2`;
   - le disuguaglianze (`vars(3) > xA; ...`) sono inserite nel vettore dei residui di `fsolve` come valori logici 0/1. Non sono vincoli: il solutore tende anzi a renderle false;
   - il punto iniziale è casuale e il ciclo termina quando lo scarto massimo è sotto **0.40 m**, una tolleranza enorme; può anche non terminare mai.
5. **Parametri macchina scritti a mano e diversi tra gli script.** Si trovano VNS 2024, VNS 2026, DEMO, PROTO, CEFTR e PILOT come blocchi commentati. `_n_wires` usa DEMO; B0 vale 5.6 negli script e 5.7 nell'Excel di MADE. Anche la sezione WP 0.46 × 0.33 m e `Iop` sono arbitrari.
6. **Geometria del pacchetto nei file emag.**
   - Filamenti su una griglia uniforme NL×NT, non i turn reali di MADE.
   - Lo scostamento toroidale è angolare (`phi' = phi0 + v/r`), mentre la bobina è piana: lo scostamento va fatto lungo la normale al piano della bobina.
   - `B_TF` è moltiplicato per un fattore empirico `*1.08`.
7. **Biot–Savart.**
   - Regola del punto medio con regolarizzazione di Hurwitz di raggio `a = dy_WP/2`, circa 165 mm (base) o `sqrt(A_cond/pi)` (full).
   - Con segmenti lunghi quanto la distanza tra conduttori vicini, la regola del punto medio non è accurata nel campo vicino.
   - Nella versione base le sorgenti sono i punti iniziali dei segmenti, mentre i punti di valutazione sono etichettati come punti medi.
8. **Induttanza.**
   - In `_full` il potenziale vettore usa i `dl` dei punti di campo al posto di quelli delle sorgenti (`sum(DLx(idx_p)./R(idx_p,:),2)`).
   - Il ciclo gira solo su `q = 1`.
   - `N_spire` e `I_spira` non sono definiti, quindi lo script si ferma con un errore.
   - La matrice `R` ha dimensione di circa 130 000 × 10 000, cioè circa 11 GB.
   - Nella versione base `Wmag` omette N² (l'energia `E` con N² è invece corretta).
9. **Altro.**
   - L'ultima cella di `TF_shape_opt_BF_Ncoils.m` legge un file con data fissa (`'TF shape - VNS 07 2026 - 09-Jul-2026.txt'`).
   - `memory` esiste solo su Windows (c'è un fallback).
   - `delete` di file che potrebbero non esistere.

### 5.3 Piano proposto

Il piano non è ancora implementato. Nuova cartella `coil3d/`, eseguita dopo la scelta della soluzione e dopo il passo 5b del main. Tutto è letto da `row` e `p`, nessun parametro è scritto a mano.

**1. Geometria dalla soluzione MADE** (`tf3d_from_design.m`, driver)

- Centroide della corrente sulla gamba interna: media dei raggi dei turn, con le stesse ricorsioni di §1 (oppure le uscite `x`, `y` di `compute_discrete_field_profile`). Distanza dalla faccia del case: `d_c = Ri_ - r_c`.
- Raggi della linea di corrente:
  - `r1 = RTFi - d_c`, sulla gamba interna;
  - `r2 = RTFo + d_c`, sulla gamba esterna, dove `RTFo = (R0 + R0/A)*(1/ripple)^(1/n_TF)` è la faccia lato plasma, e il WP è specchiato.
- Corrente reale per bobina: `Iop × sum(n_turns)`, maggiore o uguale a `g.NI` per l'arrotondamento.
- Nose della gamba esterna: `nose_ol = tf3d_nose_ol_factor × nose_il`, con default 0.5 come negli script originali. Da confermare con l'utente.

**2. Forma bending free** (`tf_bending_free_shape.m`), senza ANSYS

- Parametrizzazione con t ∈ [−π/2, π/2]:
  ```
  intB(r) = ∫_{r1}^{r} B_eff dr' = k (1 − sin t)
  k = intB(r2)/2
  dz = −k sin t / B_eff dt
  ```
  Con z = 0 in r2, più la gamba dritta in r1.
- `B_eff` è il campo toroidale medio sul pacchetto. Il vincolo è f_n = NI·B_eff, e T = (NI/2)∫B_eff dr è costante.
- Primo passo: `B_eff = μ0 n_TF NI / (4π r)`, che dà la D di Princeton (File).
- Passi successivi: `B_eff(r)` ricavato dalla forza normale reale per unità di lunghezza, calcolata con il Biot–Savart 3D dei turn MADE sulla forma corrente (bobine vicine raggruppate, vedi punto 5). Iterare finché lo spostamento relativo è sotto `tf3d_shape_tol`.
- Controllo: con n_TF grande il risultato deve tendere alla D analitica.

**3. Fit a 3 archi deterministico** (`tf_three_arc_fit.m`)

- Arco 1: tangente verticale in A (sommità della gamba dritta), centro `(xA + r1, zA)`.
- Arco 3: tangente verticale in B = (r2, 0), centro `(xB − r3, 0)`.
- Arco 2: tangente internamente a entrambi, `|C2 − C1| = R2 − r1` e `|C2 − C3| = R2 − r3`. Dati r1, r3 e R2, il centro C2 si ricava come intersezione di due circonferenze.
- Tre incognite (r1, R2, r3), minimizzando l'RMS della distanza punto-arco dalla curva bending free con `fminsearch`. Punto iniziale da una griglia grossolana.
- Uscite: centri, raggi, angoli, punti di tangenza, scarto massimo. Lo scarto va riportato e non usato come condizione di ciclo.

**4. Curve di offset ed export CAD** (`export_tf3d_shape.m`)

- Archi concentrici dalla linea di corrente. Verso l'esterno della D: faccia posteriore del WP, poi case esterno con il nose. Verso l'interno: faccia lato plasma del WP, poi case lato plasma.
- Stesso formato dell'originale, in mm: `r_case_CL z_case_CL r_case_i ...`, più un file con i dettagli degli archi.

**5. Modello 3D a filamenti** (`tf3d_coil_filaments.m`, `biot_savart_segments.m`)

- Ogni turn MADE diventa un filamento chiuso che segue la linea di corrente con offset `u` nel piano (lungo la normale uscente dalla D) e `v` lungo la **normale al piano della bobina**:
  ```
  X = (r + u n_r) e_r(φ0) + v e_φ(φ0) + (z + u n_z) e_z
  u = r_c − y_turn   (positivo lontano dal plasma)
  v = x_turn
  ```
- Bobina 1: tutti i turn. Bobine 2..n_TF: un filamento per layer, oppure 2 filamenti di Gauss a ±w/(2√3) per layer, per il costo di calcolo.
- **Biot–Savart esatto per segmenti rettilinei**:
  ```
  B = μ0 I/(4π) · (|a|+|b|) (a×b) / (|a||b| (|a||b| + a·b))
  ```
  dove a e b vanno dal punto agli estremi del segmento. Contributo nullo se il punto giace sulla retta del segmento. Serve un calcolo a blocchi per la memoria.
- Punti di valutazione: centroidi dei turn della bobina 1 lungo il percorso, solo metà superiore (simmetria z).
- **Picco 3D per turn:** `B_center_3D(s) + (B_discrete_2D − B_center_2D)`, cioè il campo 3D al centro più l'incremento locale di auto-campo del modello 2D validato.

**6. Uscite da confrontare con MADE**

| Uscita 3D | Confronto con MADE / controllo |
|---|---|
| Media toroidale di B_φ in (R0, z=0) | Deve dare **esattamente** μ0 n_TF NI/(2π R0) (legge di Ampère): test di correttezza del Biot–Savart. Confrontare con p.B0 |
| Ripple a R0 + a sul piano equatoriale, (Bmax−Bmin)/(Bmax+Bmin) | Contro `p.ripple` (MADE ricava RTFo da una formula approssimata) |
| Picco di campo 3D per turn, massimo e posizione (gamba dritta / curva / gamba esterna) | Contro `B_discrete` 2D. A metà della gamba interna 3D e 2D devono coincidere entro qualche per cento |
| Tensione T = F_z(metà superiore)/2 dalla forza di Lorentz sui filamenti | Contro `T_bf = 0.5 k_bf n_TF NI² μ0/2π`, la forza assiale usata dal modello FE 2D |
| Forza di centraggio per bobina | Informativa, per il case e il wedging |
| Energia e induttanza equivalente per bobina, L_eq = 2W/(Iop² n_TF): formula di Neumann tra le linee di corrente; per la mutua sulla stessa bobina, distanza regolarizzata con la GMD della sezione rettangolare, circa 0.2235(w+h) | Contro l'induttanza a guscio della scansione (`L = μ0 (Ntot k_bf)² r_bf/2 (I0+2I1+I2)/n_TF`), che determina τ di scarica e quindi l'hot spot |

**7. Grafici** (`postprocess/plot_tf3d.m`)

- forma: bending free, 3 archi, offset di case e WP;
- vista 3D con |B| sui turn;
- picco per turn lungo l'ascissa curvilinea, con il riferimento 2D;
- B_φ(R) sul piano equatoriale con il ripple.

**8. Parametri Excel** (nuova categoria "3D coil")

| Parametro | Default | Significato |
|---|---|---|
| `tf3d_shape_mode` | 1 | 0 = D analitica, 1 = iterata con il campo 3D |
| `tf3d_n_seg` | 360 | Segmenti per bobina |
| `tf3d_shape_tol` | 1e-4 | Tolleranza dell'iterazione di forma |
| `tf3d_max_iter` | 20 | Iterazioni massime |
| `tf3d_nose_ol_factor` | 0.5 | Nose gamba esterna / nose gamba interna |
| `tf3d_use_arcs` | 1 | 1 = modello 3D costruito sulla forma a 3 archi |
| `tf3d_export` | 1 | Scrive i file per CAD |

**9. Validazione** (`validation/validate_tf3d.m`)

- spira circolare, campo sull'asse: B = μ0 I R² / (2(R²+z²)^{3/2});
- controllo di Ampère a R0;
- limite n_TF grande della forma iterata contro la D analitica;
- 3D contro 2D a metà della gamba interna sul design 7. La geometria del design 7 è in `validation/validate_mech_surrogate_2026.m`, caso B.

**10. Main:** passo 6 facoltativo "Run the 3D coil step (shape + 3D field)? [y/N]"; l'export ANSYS diventa il passo 7.

**11. Documentazione:** aggiungere una sezione al report e al file findings, con l'elenco dei problemi degli originali corretti (§5.2).

### 5.4 Domande da fare all'utente

1. Nel 3D il WP va costruito sulla forma a 3 archi (quella costruita davvero) o sulla forma bending free?
2. Confermare `nose_ol = 0.5 × nose_il` e che la piastra lato plasma della gamba esterna sia uguale a quella della gamba interna.
3. Servono anche i file per ANSYS 3D (per esempio i punti della linea di corrente per SOURC36 o per un modello solido), oltre all'export CAD?

---

## 6. File chiave e dimensioni

| File | Righe circa | Ruolo |
|---|---|---|
| `main_WP_TF_design.m` | 110 | Punto d'ingresso |
| `search/scan_wp_designs.m` | 440 | Scansione |
| `physics/wp_mech_surrogate.m` | 1400 | FE 2D |
| `physics/wp_field_at_points.m` | 90 | Campo 2D |
| `physics/compute_discrete_field_profile.m` | 150 | Picco per turn |
| `validation/validate_mech_surrogate_2026.m` | 110 | Validazione PASS/FAIL |
| `docs/REPORT_codice_e_modello_meccanico.md` | 415 | Guida completa |

I file nella radice con nomi come `WP_TF_*_Design_*.m`, `PLOT_TF_*.m`, `CICC_*.m` sono **legacy**: non sono usati dalla pipeline e non vanno modificati.
