# Modello analitico a strati della gamba interna (jacket e case accoppiati)

Diario del lavoro: equazioni, ipotesi, codice, confronti con il FE 2D e
decisioni. Codice: `physics/wp_layered_cylinder.m`; confronto:
`validation/verify_layered_cylinder.m`; dati di riferimento:
`validation/results/mech_reference_designs.mat` (13 design),
`validation/results/jacket_surrogate_calibration.csv` (jacket, 23 run FE),
`validation/results/mech_reference_fe.csv` (jacket e case, run FE dedicati,
`validation/run_mech_reference_fe.m`).

## 1. Perché

Oggi jacket e case sono dimensionati da due modelli separati (surrogato per
il jacket, formula del cilindro spesso per il vault) con criteri diversi, e
si scambiano solo S_z e l'altezza del WP
(`docs/DIMENSIONAMENTO_JACKET.md`). L'obiettivo è un modello fisico unico,
veloce come una formula, in cui:
- il carico radiale si accumula verso il nose (nel FE +70% sul jacket dal
  lato plasma ai layer profondi);
- il carico circonferenziale del wedging si ripartisce tra nose, pareti
  laterali e WP secondo le rigidezze (un jacket più spesso scarica il nose);
- jacket e case si verificano con lo stesso criterio (Pm ≤ Sm, Pm+Pb ≤ 1.5 Sm).

## 2. Ipotesi

1. **Cilindro a 360°.** Le bobine sono incuneate: i fianchi dei case si
   toccano e trasmettono la compressione circonferenziale. Per il carico
   centripeto la gamba interna lavora come un cilindro continuo; ogni settore
   (2π/n_TF) è rappresentato dalla sua corda 2 r tan(π/n_TF).
2. **Simmetria assiale, deformazione piana generalizzata.** Incognite: lo
   spostamento radiale u(r) e una deformazione assiale uniforme ε₀.
3. **Strati concentrici omogenei**, dal foro verso l'esterno: nose (acciaio, dal foro ad arco R_b a R_j, §3.0), isolante di terra, i layer del WP separati dall'isolante
   tra layer, isolante di terra, piastra lato plasma (acciaio, spessore
   dr_plasma_side). Ogni anello di WP contiene, lungo la corda del settore,
   le due pareti laterali del case, i due isolanti di terra laterali e i turn
   del layer.
4. **Caso primario**: forza di Lorentz e forza verticale T_bf, senza
   raffreddamento (come il criterio Pm, Pm+Pb).

## 3. Geometria del settore e rigidezze equivalenti

### 3.0 Geometria come nel FE (e in ANSYS GPS)

- **Foro ad arco.** Il FE (come il modello ANSYS, `WEDGE_INTR_INS`) usa per
  il foro l'arco di raggio R_b = R_k / cos(π/n_TF), che passa per gli
  spigoli del foro piano a distanza R_k. Al centro il nose è quindi spesso
  R_j − R_b e non R_j − R_k: con n_TF = 12 sono 17–25 mm in meno, il
  15–35% del nose (tabella in §8.3). Default `bore = 'arc'`; `'flat'` usa R_k.
- **Isolante di wedge** (t_w = 1.5 mm, `wedge_insulation`) su entrambi i
  fianchi: è in serie nel percorso circonferenziale di ogni anello, nose e
  piastra compresi. Per un anello di acciaio la cedevolezza circonferenziale
  diventa S₂₂ = 1/E + f_w/E_ins, con f_w = 2 t_w / corda.
- **Cavità a inviluppo convesso**: la cavità del FE è l'inviluppo convesso
  dei layer allargato dell'isolante di terra. Dove un layer è più stretto di
  quello sopra, lo scalino è riempito di isolante e non di acciaio del case.
  Default `cavity = 'hull'`; `'exact'` usa la larghezza di ogni layer.

### 3.1 Cella del turn (isolante avvolgente)

Cella Cond_w × Cond_h: il conduttore incamiciato (a × b, con
a = Cond_w − 2 t_ins e b = Cond_h − 2 t_ins) è avvolto dall'isolante di turn
t_ins. Dentro c'è il cavo SC_w × SC_h con il jacket di spessore JT. L'isolante
ha due moduli: E_ins attraverso lo spessore (12 GPa) ed E_ins,t nel piano
(20 GPa). r_SC e il raccordo esterno non entrano nelle rigidezze.

Conduttore incamiciato:

    radiale:          E_mid,r = 1 / ( 2 JT/(b E_j) + SC_h/(b E_cbl) )     colonna centrale: pareti sup./inf. e cavo in serie
                      E_jc,r  = ( 2 JT E_j + SC_w E_mid,r ) / a              pareti laterali in parallelo
    circonferenziale: E_mid,t = 1 / ( 2 JT/(a E_j) + SC_w/(a E_cbl) )
                      E_jc,t  = ( 2 JT E_j + SC_h E_mid,t ) / b

Turn: l'isolante trasversale al carico è in serie con tutto il conduttore,
quello parallelo al carico è in parallelo:

    radiale:          E_col,r = 1 / ( 2 t_ins/(Ch E_ins) + b/(Ch E_jc,r) )
                      E_r,turn = ( 2 t_ins E_ins,t + a E_col,r ) / Cw
    circonferenziale: E_band,t = 1 / ( 2 t_ins/(Cw E_ins) + a/(Cw E_jc,t) )
                      E_t,turn = ( 2 t_ins E_ins,t + b E_band,t ) / Ch
    assiale:          E_z,turn = ( A_j E_j + A_cbl E_cbl + A_ins E_ins,t ) / A_cella

La prima versione (`cell = 'open'`) metteva le pareti del jacket in
parallelo su tutta la cella, scavalcando l'isolante di turn. Così
sovrastimava la rigidezza circonferenziale del WP, che sottraeva carico al
nose e alla piastra (§8).

Nota: `cavo_stiffness` (in `size_cicc_cable`) usa ancora un'altra
combinazione, che mescola le frazioni di larghezza e di altezza della cella.

### 3.2 Anello

Con f_k la frazione della corda occupata dal componente k:
- pareti laterali del case: w_wall = corda − W_cav − 2 t_w;
- isolante: isolante di terra laterale, riempimento dell'inviluppo, wedge;
- WP: W = n_turn Cond_w.

    radiale        E_r = Σ f_k E_r,k            (componenti in parallelo)
    circonferenz.  1/E_t = Σ f_k / E_t,k         (componenti in serie)
    assiale        E_z = Σ f_k E_z,k             (parallelo)

L'isolante laterale ha lo spessore in direzione toroidale: è caricato
attraverso lo spessore (E_ins) nella direzione circonferenziale e nel piano
(E_ins,t) in quella radiale e assiale. Gli anelli di solo isolante (di terra,
tra layer) hanno E_r = E_ins ed E_θ = E_z = E_ins,t.

Gli anelli di acciaio (nose, piastra) sono isotropi con l'accoppiamento di
Poisson (ν = 0.29), più l'isolante di wedge in serie; gli anelli misti
sono ortotropi senza accoppiamento di Poisson. Rigidezze in deformazione
piana generalizzata, ε = (ε_r, ε_θ, ε_z):

    σ = C ε,   acciaio: C = S⁻¹, S = (1/E) [1 −ν −ν; −ν 1 −ν; −ν −ν 1] + diag(0, f_w/E_ins, 0)
               misto:   C = diag(E_r, E_t, E_z)

I moduli si leggono da p in GPa, come in `wp_mech_surrogate`.

## 4. Equazioni in un anello

Congruenza: ε_r = du/dr, ε_θ = u/r, ε_z = ε₀.

Equilibrio radiale con la forza di Lorentz per unità di volume f(r):

    dσ_r/dr + (σ_r − σ_θ)/r + f(r) = 0

Forza di Lorentz (smeared). In un anello di WP la densità di corrente è
j = n_TF n_k I_op / (π (r₂² − r₁²)); il campo, per la simmetria assiale, è
B_θ(r) = μ₀ I_enc(r) / (2π r), con I_enc la corrente racchiusa entro r.
Nell'anello: I_enc(r) = I_0 + j π (r² − r₁²), con I_0 la corrente degli anelli
più interni. La forza è centripeta:

    f(r) = −j B_θ(r) = a₁ r + a₋₁ / r,   a₁ = −μ₀ j²/2,   a₋₁ = −μ₀ j (I_0 − j π r₁²)/(2π)

Sostituendo σ = C ε nell'equilibrio (i termini in C₁₂ si cancellano):

    C₁₁ (u'' + u'/r) − C₂₂ u / r² = −f(r) − (C₁₃ − C₂₃) ε₀ / r

Soluzione:

    u(r) = A r^k + B r^(−k) + P r³ + Q r + ε₀ Q_e r,     k = √(C₂₂/C₁₁)
    P = −a₁ / (9 C₁₁ − C₂₂)
    Q = −a₋₁ / (C₁₁ − C₂₂),   Q_e = −(C₁₃ − C₂₃)/(C₁₁ − C₂₂)

Se C₁₁ = C₂₂ (anello isotropo, k = 1) il termine in r diventa r ln r con
coefficienti −a₋₁/(2 C₁₁) e −(C₁₃ − C₂₃)/(2 C₁₁); per un anello isotropo
C₁₃ = C₂₃ e il termine in ε₀ sparisce (Lamé: u = A r + B/r + P r³ + …).

## 5. Collegamento degli anelli

Incognite: A_i, B_i per ogni anello (2N) e ε₀. Equazioni:

- σ_r = 0 sul foro del nose (r = R_k) e sulla faccia lato plasma (r = R_i);
- a ogni interfaccia: u continuo e σ_r continuo (2(N−1) equazioni);
- equilibrio assiale: Σ_i ∫ σ_z 2π r dr = n_TF T_bf (integrale con Gauss a
  6 punti per anello).

Sono 2N + 1 equazioni lineari; con 13–28 layer, nose, piastra e isolanti
N ≈ 30–60. La soluzione è esatta per il modello (residuo ~10⁻¹⁶) e
richiede meno di un millisecondo. Controllo indipendente: con le superfici
libere ∫σ_θ dr = ∫ f r dr su tutta la sezione (rispettato a 10⁻⁴, errore
della quadratura a trapezi).

## 6. Dall'anello alle pareti del jacket e del case

Per il layer k, a metà altezza, dalla soluzione si hanno ε_r, ε_θ, σ_θ
dell'anello e ε₀.

1. **Circonferenziale.** La tensione circonferenziale è la stessa in tutti
   i componenti in serie dell'anello. Quindi:
   - deformazione del turn: ε_θ,turn = σ_θ / E_t,turn;
   - tensione della fascia: σ_band = E_band,t ε_θ,turn;
   - deformazione del conduttore incamiciato: ε_jc,θ = σ_band / E_jc,t.
2. **Radiale.** Deformazione del turn: ε_r = −p_k / E_r,turn (colonna, §7)
   oppure ε_r dell'anello. Poi:
   - tensione della colonna: σ_col = E_col,r ε_r;
   - deformazione del conduttore incamiciato: ε_jc,r = σ_col / E_jc,r.
3. **Pareti.** Ogni parete del jacket ha una deformazione imposta nella
   direzione in cui è "in parallelo", una tensione imposta nella direzione
   in cui è "in serie", e ε_z = ε₀. Le altre due tensioni vengono dalla
   legge di Hooke 3D:

        ε_par = (σ_par − ν (σ_ser + σ_z))/E,   ε₀ = (σ_z − ν (σ_par + σ_ser))/E

   - pareti laterali: ε_par = ε_jc,r (radiale), σ_ser = E_mid,t ε_jc,θ;
   - pareti superiore e inferiore: ε_par = ε_jc,θ, σ_ser = E_mid,r ε_jc,r;
   - pareti laterali del case: ε_par = ε_r dell'anello, σ_ser = σ_θ dell'anello.

4. **Verifica.**
   - Pm del jacket del layer = massima Tresca tra parete laterale e parete
     superiore.
   - Il nose e la piastra si linearizzano attraverso lo spessore come le
     linee di classificazione del FE: membrana = media, flessione = parte
     lineare di σ_r, σ_θ, σ_z.

L'effetto Poisson su σ_z è essenziale: con σ_z = E ε₀ il jacket risultava
~100 MPa troppo caricato assialmente (289 contro ~190 MPa sul design 10).

## 7. Varianti del percorso radiale del WP

- **`radial` = 'ring'**: il WP è incollato al case, il carico radiale si
  ripartisce con le pareti laterali (rigide) nell'anello.
- **`radial` = 'column'** (default): il WP è una colonna che scorre sulle
  pareti (nel FE l'interfaccia WP/case è un contatto con attrito); la
  pressione radiale sul layer k è la forza di Lorentz dei layer sopra per
  unità di larghezza di ogni layer (dove un layer è più stretto, i turn in
  sbalzo appoggiano sulla spalla del case):

        p_k = Σ_{j<k} F_j/W_j + F_k/(2 W_k)

  con attrito sulle pareti (Janssen): la compressione circonferenziale σ_θ
  preme il WP sulle due pareti e ogni layer di altezza h_k cede al case fino a
  2 μ |σ_θ| h_k per unità di lunghezza.
  - Ricorsione: p_in = 0 in cima e per ogni layer
    Δp = F_k/W_k − 2 μ |σ_θ,k| h_k / W_k,
    p_k = max(p_in + Δp/2, 0), poi p_in ← max(p_in + Δp, 0).
  - μ = `mu_case` è un coefficiente **efficace**, l'unico parametro tarato
    sul FE: vale 0.1, mentre il contatto del FE ha 0.2. La pressione sulle
    pareti non è ovunque |σ_θ| (raccordi, scalini, layer stretti).
  - Deformazione radiale del turn: ε_r = −p_k / E_r,turn.
  - La colonna cambia solo il recupero del jacket. Nose, piastra e pareti del
    case vengono dalla soluzione ad anelli.
- **`wall_radial`**: quota delle pareti laterali nella rigidezza radiale
  dell'anello (1 incollato, 0 scorrevole).

## 8. Risultati (02/10/2026)

### 8.1 Jacket, 23 run FE (Pm primario per layer)

Le prime righe sono la storia dello sviluppo. L'ultima riga del modello è la
versione attuale (default).

| Versione | Rapporto medio modello/FE | rms | 90° perc. di \|r−1\| | Massimo del design, min … max |
|---|---:|---:|---:|---|
| anello, cella 'open', foro piano, σ_z = E ε₀ | — | — | — | 15–20% sopra il FE |
| anello, cella 'open', foro piano, Poisson | 0.86 | 18.3% | — | 0.62 … 1.08 |
| colonna μ = 0, cella 'open', foro piano | 1.04 | 14.1% | — | 0.88 … 1.29 |
| anello, geometria FE, cella 'wrapped' | 0.88 | 18.8% | — | 0.59 … 0.90 |
| colonna μ = 0.2, geometria FE, 'wrapped' | 0.92 | 13.1%* | 0.18* | 0.71 … 0.93* |
| colonna μ = 0.05 | 1.02 | 9.3%* | 0.15* | 0.87 … 1.18* |
| **colonna μ = 0.1 (default)** | **0.98** | **8.3%*** | **0.135*** | **0.85 … 1.09*** (rms 9.3%) |
| colonna μ = 0.15 | 0.95 | 10.3%* | 0.17* | 0.74 … 1.01* |
| surrogato (riferimento) | 1.0 | 7% | — | LOO 0.87 … 1.23 |

\* senza le due run FE anomale escluse anche dal fit del surrogato
(s1_1@jt+1.0, s2_4@jt+1.0). Con tutte le 430 righe il default dà: medio
0.983, rms 8.8%, massimo del design 0.70 … 1.09.

**Risposta allo spessore del jacket** (Pm massimo):

| Design | FE +0.5 mm / +1.0 mm | prima versione | default attuale |
|---|---|---|---|
| bench | −4.4% / −8.4% | −4.8% / −9.3% | −7.3% / −11.8% |
| d7 | −6.4% / −11.0% | — | −9.0% / −16.2% |
| d10 | −5.5% / −10.9% | −6.9% / −12.9% | −10.7% / −18.8% |

La versione attuale è circa 1.5 volte più sensibile del FE allo spessore del
jacket. Nel dimensionamento sarebbe ottimista quando si ingrossa il jacket.
È un limite da tenere presente prima di usarla nel loop.

### 8.2 Componenti delle pareti (design 7, turn centrali, prima versione)

| Layer | FE parete lat. σ_r / σ_θ / σ_z | modello (anello) | FE parete sup. σ_θ / σ_r / σ_z | modello (anello) |
|---|---|---|---|---|
| L1 | +67 / −106 / 277 | +67 / −24 / 256 | −283 / 1 / 207 | −382 / 0 / 133 |
| L7 | −260 / −187 / 159 | −28 / −24 / 229 | −315 / −132 / 159 | −378 / −6 / 132 |
| L13 | −289 / −234 / 137 | −97 / −23 / 209 | −398 / −159 / 127 | −372 / −10 / 133 |

FE = caso totale, con raffreddamento; MPa. La cella in serie e parallelo
separa troppo i percorsi: nel FE le pareti laterali portano anche
compressione circonferenziale e le pareti superiori anche compressione
radiale. Il trasferimento passa per i raccordi arrotondati e per il cavo
confinato. La colonna con attrito corregge il percorso radiale in media, non
la ripartizione locale.

### 8.3 Case (Pm primario, MPa; FE / modello)

| Versione | Design | nose al centro | piastra lato plasma | parete laterale 50% |
|---|---|---|---|---|
| prima (foro piano, cella 'open') | bench | 655 / 521 (−20%) | 545 / 437 (−20%) | 402 / 341 |
| | d7 | 606 / 506 (−17%) | 534 / 407 (−24%) | 357 / 305 |
| + foro ad arco, wedge, inviluppo | bench | 655 / 586 (−11%) | 545 / 500 (−8%) | 402 / 349 |
| | d7 | 606 / 533 (−12%) | 534 / 437 (−18%) | 357 / 328 |
| + cella 'wrapped' (default) | bench | 655 / 601 (−8%) | 545 / 514 (−6%) | 402 / 366 (−9%) |
| | d7 | 606 / 553 (−9%) | 534 / 452 (−15%) | 357 / 324 (−9%) |

**Tutte le 19 run FE di riferimento** (versione di default;
`validation/results/mech_reference_fe.csv`, `run_mech_reference_fe`, tutte
valide). Spessore del nose al centro come nel FE (mm), poi FE / modello in MPa.

| Run | nose FE [mm] | nose al centro | piastra lato plasma | parete laterale 50% |
|---|---:|---|---|---|
| bench | 111 | 655 / 601 (−8%) | 545 / 514 (−6%) | 402 / 366 (−9%) |
| d7 | 107 | 606 / 553 (−9%) | 534 / 452 (−15%) | 357 / 324 (−9%) |
| d10 | 154 | 631 / 599 (−5%) | 526 / 522 (−1%) | 403 / 436 (+8%) |
| s1_1 | 114 | 577 / 505 (−12%) | 448 / 390 (−13%) | 321 / 269 (−16%) |
| s1_13 | 119 | 566 / 497 (−12%) | 410 / 375 (−9%) | 312 / 264 (−15%) |
| s2_4 | 35 | 782 / 583 (−25%) | 463 / 365 (−21%) | 318 / 248 (−22%) |
| s1_2 | 107 | 589 / 514 (−13%) | 463 / 396 (−14%) | 326 / 272 (−17%) |
| s1_4 | 119 | 561 / 494 (−12%) | 412 / 373 (−9%) | 309 / 250 (−19%) |
| s1_7 | 75 | 649 / 538 (−17%) | 468 / 396 (−15%) | 332 / 273 (−18%) |
| s1_9 | 101 | 585 / 507 (−13%) | 425 / 376 (−12%) | 319 / 253 (−21%) |
| s1_17 | 109 | 576 / 496 (−14%) | 396 / 370 (−7%) | 314 / 261 (−17%) |
| s2_1 | 33 | 788 / 573 (−27%) | 491 / 373 (−24%) | 324 / 262 (−19%) |
| s2_2 | 34 | 829 / 602 (−27%) | 492 / 385 (−22%) | 327 / 259 (−21%) |
| bench@jt+0.5 | 111 | 640 / 580 (−9%) | 526 / 493 (−6%) | 392 / 352 (−10%) |
| bench@jt+1.0 | 111 | 626 / 561 (−10%) | 507 / 473 (−7%) | 383 / 339 (−11%) |
| d7@jt+0.5 | 107 | 590 / 530 (−10%) | 511 / 429 (−16%) | 347 / 311 (−10%) |
| d7@jt+1.0 | 107 | 576 / 511 (−11%) | 490 / 409 (−17%) | 337 / 298 (−12%) |
| d10@jt+0.5 | 154 | 613 / 574 (−6%) | 504 / 496 (−2%) | 391 / 415 (+6%) |
| d10@jt+1.0 | 154 | 597 / 552 (−8%) | 483 / 473 (−2%) | 382 / 397 (+4%) |

Rapporto modello/FE: nose 0.73 … 0.95 (medio 0.87), piastra 0.76 … 0.99,
parete laterale 0.78 … 1.08.

Osservazioni:
- **L'errore sul nose cresce quando il nose si assottiglia**:
  - −5% con 154 mm (d10);
  - −8 … −14% tra 100 e 120 mm;
  - −17% con 75 mm (s1_7);
  - −25 … −27% con 33–35 mm (s2_x).

  Non è quindi un fattore costante. Ipotesi da verificare sulle componenti
  (σ_θ e σ_z separati nel nose del FE, che le run non salvano):
  - la flessione del vault come trave tra le pareti, che nel FE dà anche
    membrana sulla linea centrale;
  - una ripartizione assiale diversa nel nose sottile.
- **Spessore del jacket e case.** Sul bench, +0.5 / +1.0 mm di jacket danno:
  - nose: FE −2.3% / −4.4%, modello −3.5% / −6.7%;
  - piastra: FE −3.5% / −7.0%, modello −4.1% / −8.0%.

  Il FE conferma che un jacket più spesso scarica il case, ma poco (~2% ogni
  0.5 mm sul nose). Il modello coglie il segno ed esagera l'entità, come per
  il jacket.
- **Pareti laterali** sottostimate del 15–22% sui design della scansione
  tarata. Il modello di parete (deformazione radiale dell'anello, σ_θ
  dell'anello) è il più grezzo.

Spessore del nose al centro: R_j − R_k contro il FE R_j − R_k/cos(π/n_TF).

| Design | R_j − R_k | FE | Design | R_j − R_k | FE |
|---|---:|---:|---|---:|---:|
| bench | 135.4 | 110.7 | s1_7 | 95.0 | 74.7 |
| d7 | 133.0 | 107.2 | s1_9 | 119.0 | 101.4 |
| d10 | 179.0 | 153.7 | s1_17 | 124.0 | 109.4 |
| s1_1 | 135.0 | 114.4 | s2_1 | 51.0 | 32.9 |
| s1_13 | 136.0 | 119.2 | s2_2 | 51.0 | 34.0 |
| s2_4 | 51.0 | 35.1 | s1_2 | 128.0 | 106.7 |
| s1_4 | 137.0 | 119.0 | | | |

(mm)

**Conseguenza per lo scan.** `size_case_vault` dimensiona il nose come
R_j − R_k, cioè con il foro piano. Se il case reale (o il modello ANSYS) ha il
foro ad arco, il nose al centro è più sottile di 17–25 mm, e i vault
dimensionati dallo scan risultano meno verificati di quanto lo scan creda.
Va deciso quale geometria è quella di progetto. Se è l'arco, R_k va
interpretato come raggio dello spigolo, oppure il vault va verificato su
R_j − R_k/cos(π/n_TF).

## 9. Stato e prossimi passi

- **Jacket.** Con un solo parametro tarato (μ efficace = 0.1) il modello è
  vicino al surrogato: rms 8.3% contro 7%, massimo del design 0.85 … 1.09.
  Però sovrastima il beneficio di un jacket più spesso (§8.1).
- **Case** (19 run FE). Il nose è sottostimato del 5–27% (medio −13%), di più
  quando è sottile; la piastra dell'1–24%; la parete laterale tra −22% e +8%.
  Il criterio concordato (±5%) non è rispettato. Prima di correggere serve il
  confronto per componenti nel nose sottile (s2_4).
- **Prossimo passo.** Una cella elementare risolta numericamente una volta
  per grade: un turn con raccordi, cavo, jacket e isolante, con le tre
  deformazioni medie unitarie. Dà le rigidezze equivalenti esatte e le
  matrici di localizzazione delle pareti, e serve per la sensibilità a JT e
  per il trasferimento tra pareti laterali e superiori (§8.2). L'anello
  resta analitico.
- **Decisione.** Si adotta il modello nello scan se il jacket ha rms ≤ 10%
  e 90° percentile ≤ 1.10, e il case sta entro ±5%. Altrimenti si resta con
  il surrogato del jacket e si aggiunge un surrogato del case.

## 10. Verifica affiancata con il FE (03/10/2026)

`validation/compare_layered_vs_fe.m` mette FE 2D e modello sullo stesso
design e sullo stesso carico primario, senza raffreddamento (nel FE
T_op = T_ref). Non confronta solo le tensioni di verifica ma **dove va il
carico**. Per questo `wp_mech_surrogate` restituisce ora anche, per ogni
elemento, centroide, area, materiale e tensione media (`out.elem`); la
soluzione FE non cambia.

### 10.1 Forza circonferenziale attraverso il piano x = 0 (nose / WP / piastra)

| Design (nose) | FE | modello | totale FE / modello [MN/m] |
|---|---|---|---|
| d10 (154 mm) | 69 / 26 / 5 % | 67 / 27 / 6 % | 93.0 / 92.4 |
| bench (111 mm) | 60 / 33 / 7 % | 56 / 37 / 7 % | 85.4 / 85.5 |
| s1_7 (75 mm) | 41 / 54 / 5 % | 34 / 61 / 5 % | 103.1 / 102.2 |
| s2_4 (35 mm) | 25 / 70 / 6 % | 18 / 77 / 5 % | 109.2 / 109.1 |

### 10.2 Membrana nel case (MPa, FE / modello)

| Design | nose circ. | nose ass. | piastra circ. | piastra ass. | parete radiale | parete ass. |
|---|---|---|---|---|---|---|
| d10 | −416 / −405 | 216 / 194 | −249 / −281 | 277 / 241 | +43 / −191 | 349 / 246 |
| bench | −464 / −432 | 191 / 169 | −292 / −296 | 253 / 219 | +62 / −130 | 339 / 236 |
| s1_7 | −558 / −464 | 92 / 74 | −285 / −253 | 183 / 144 | +29 / −72 | 252 / 164 |
| s2_4 | −774 / −568 | 7 / 15 | −313 / −252 | 150 / 113 | +38 / −70 | 227 / 132 |

### 10.3 WP, piano di mezzeria, media per layer (bench, MPa, FE / modello)

| Layer | circonferenziale | radiale | Pm jacket |
|---|---|---|---|
| L1 | −62 / −101 | −8 / −8 | 516 / 533 |
| L3 | −75 / −102 | −86 / −36 | 557 / 536 |
| L6 | −83 / −102 | −155 / −67 | 596 / 535 |
| L9 | −106 / −95 | −154 / −85 | 565 / 540 |
| L11 | −117 / −123 | −139 / −87 | 594 / 550 |

### 10.4 Lettura

1. **Effetto telaio assente.** Nel FE le pareti laterali sono tese
   radialmente (+30…+60 MPa): il WP spinge sul nose, il nose cede e le
   pareti lo appendono alla piastra. Nel modello le pareti condividono la
   deformazione radiale del WP e sono compresse (−70…−190 MPa).
   - La loro tensione assiale nel FE è più alta del ~40%.
   - È un limite dell'ipotesi assialsimmetrica, non un parametro da tarare.
2. **WP troppo rigido in direzione circonferenziale nei layer esterni**
   (−100 contro −60…−80 MPa). Così il WP sottrae forza circonferenziale al
   nose, di più quando il WP è grande e il nose sottile.
   - La variante `hoop = 'centre'` (rigidezza circonferenziale del solo WP)
     corregge la forma del profilo ma sposta poco la ripartizione: il nose
     di s2_4 passa da −568 a −594 contro −774 MPa.
   - Resta da migliorare la cella del turn (raccordi, contatto solo sui
     tratti piatti).
3. **Pressione radiale al centro del WP circa il doppio della colonna del
   modello** (−155 contro −67 MPa).
   - La colonna è una media sulla larghezza: nel FE il carico si concentra
     al centro e si scarica verso le pareti.
   - Nei layer profondi la pressione cala verso il nose.
4. **Il Pm del jacket torna** (±10%) per compensazione: troppa compressione
   circonferenziale e poca radiale.

### 10.5 Proposta

- **A.** Modello analitico con effetto telaio: WP come colonna sul nose,
  case come telaio chiuso nose–pareti–piastra con il nose ad arco. È
  ricerca, con esito incerto.
- **B (consigliata).** Surrogato del jacket (§5b) e surrogato del case
  tarato sulle 19 run FE di riferimento, con il Pm del nose in funzione di
  spessore del nose, forza centripeta e larghezza del WP. Il modello a
  strati resta uno strumento di comprensione.
- **C.** Verifica del case solo col FE sulle prime soluzioni della scansione.

## 11. Scelta B: surrogato del case nella scansione (03/10/2026)

Decisione: nella scansione si usano il surrogato del jacket (§5b del
manuale) e un surrogato del case tarato sul FE. Il modello a telaio
(scelta A) resta uno studio per gli sviluppi futuri: `docs/MODELLO_TELAIO.md`.

### 11.1 Foro ad arco anche nello scan

`size_case_vault` dimensionava il nose come DTF = R_j − R_k, cioè con il
foro piano. FE e ANSYS hanno invece il foro ad arco di raggio
R_k/cos(π/n_TF). Ora DTF è lo spessore del nose al centro, sull'arco:

    R_bore = R_j − DTF,     R_k = R_bore cos(π/n_TF)

- R_k resta il parametro ANSYS (`WEDGE_INTR_INS` = R_k/cos(π/n_TF)), con lo
  stesso significato di prima.
- Nuova colonna `R_bore`. `Nose` = R_j − R_bore, ingombro radiale
  R_i − R_bore.
- L'area del case è quella del settore con il foro ad arco:
  R_i² tan(θ/2) − (θ/2) R_bore².

### 11.2 Il surrogato

    Pm,case    = (a_m + b_m ln x) · Pm,nose,LC
    Pm+Pb,case = (a_b + b_b ln x) · Pm,nose,LC,      x = h_WP / t_nose

- **Pm,nose,LC** è l'intensità di membrana del nose del modello a strati
  (opzioni di default). Porta già carichi, rigidezze e l'accoppiamento con
  il jacket.
- **x** è il rapporto tra l'altezza del WP (dalla piastra al nose) e lo
  spessore del nose al centro. Misura quanto pesa l'effetto telaio che il
  modello assialsimmetrico non vede: il rapporto FE/modello sale con
  regolarità da 1.05 (x = 2.2) a 1.37 (x = 21).
- **Bersaglio**: il massimo sulle SCL del case del FE (nose al centro, a
  metà, diagonale del vault, pareti, piastra), carico primario.

| | a | b | leave-one-design-out, modello/FE | rms |
|---|---:|---:|---|---:|
| Pm | 0.95802 | 0.13123 | 0.967 … 1.040 | 2.1% |
| Pm+Pb | 0.99183 | 0.23057 | 0.864 … 1.121 | 7.1% |

- **Dati**: le 19 run FE di riferimento (`mech_reference_rows.mat`,
  `mech_reference_fe.csv`), con x da 2.2 a 21. La validità è limitata a
  n_TF = 12 e ai cavi rettangolari; `extrapolated` segnala x fuori
  intervallo.
- **Altre forme provate** per Pm+Pb: con la larghezza del WP come secondo
  predittore il rms in leave-one-design-out non migliora (7.7%), e con
  W/corda del foro migliora poco (5.9%). Si tiene un solo predittore con un
  margine dedicato.
- **Riproducibilità**: `validation/run_case_surrogate_features.m` e
  `validation/tools/fit_case_surrogate.py`. I coefficienti valgono per le
  opzioni di default di `wp_layered_cylinder`: se cambiano, va rifatto il
  fit.

### 11.3 Nella scansione

- `case_model` = 1 (default): `size_case_surrogate` trova il DTF minimo
  sulla griglia DTF_initial + k·DTF_step, per bisezione (circa 10
  valutazioni), con:
  - case_margin·Pm ≤ Sm_case (case_margin = 1.05);
  - case_margin_PmPb·(Pm+Pb) ≤ 1.5 Sm_case (case_margin_PmPb = 1.10).
- S_z segue la nuova area del case.
- Nuove colonne: S_T_VT = Pm del case, `Case_PmPb`, `case_x`.
- `case_model` = 0: formula del vault di prima, ma con il foro ad arco.

**Accoppiamento jacket–case.** Sui design di riferimento, +0.5 mm di
jacket su tutti i layer fanno risparmiare 5 mm di nose (bench: DTF da 119 a
114 mm; d7: da 103 a 95 mm). Il surrogato lo vede perché la base è il
modello a strati.

**Scan ridotto (macchina del design 10, 22 layout da 18 turn).**

| | prima (formula del vault, foro ad arco) | surrogato del case |
|---|---|---|
| ingombro radiale della #1 | 540.1 mm | 520.9 mm |
| nose della #1 | 173 mm | 152 mm |
| Pm del case della #1 (surrogato) | — | 635 MPa |
| tempo dello scan in Octave | 169 s | 210 s |

**Verifica col FE 2D delle prime due soluzioni dello scan ridotto**
(carico primario, MPa).

| | case Pm FE / surrogato | case Pm+Pb FE / surrogato | jacket Pm FE / surrogato | jacket Pm+Pb FE / surrogato |
|---|---|---|---|---|
| #1 | 629 / 635 (+1%) | 689 / 704 (+2%) | 617 / 557 | 1028 / 900 |
| #2 | 627 / 633 (+1%) | 685 / 705 (+3%) | 614 / 553 | 1013 / 895 |

- **Case**: verificato con il margine previsto (Pm ≤ 667 MPa).
- **Jacket**: Pm+Pb supera 1.5 Sm = 1000 MPa dell'1–3%. È il layout del
  design 10, il peggiore per il surrogato del jacket (§5b del manuale): va
  corretto con `refine_jacket_fe`, come previsto dal flusso del main.
