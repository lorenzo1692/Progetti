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
3. **Strati concentrici omogenei**, dal foro verso l'esterno: nose (acciaio,
   da R_k a R_j), isolante di terra, i layer del WP separati dall'isolante
   tra layer, isolante di terra, piastra lato plasma (acciaio, spessore
   dr_plasma_side). Ogni anello di WP contiene, lungo la corda del settore,
   le due pareti laterali del case, i due isolanti di terra laterali e i turn
   del layer.
4. **Caso primario**: forza di Lorentz e forza verticale T_bf, senza
   raffreddamento (come il criterio Pm, Pm+Pb).

## 3. Rigidezze equivalenti

### 3.1 Cella del turn

Cella Cond_w × Cond_h: cavo SC_w × SC_h al centro, jacket di spessore JT,
isolante t_ins (r_SC e il raccordo esterno non entrano nelle rigidezze).

- **Radiale** (carico lungo l'altezza della cella): le pareti laterali del
  jacket e dell'isolante sono in parallelo alla colonna centrale, nella
  quale isolante, pareti superiore e inferiore del jacket e cavo sono in
  serie:

      E_mid,r = 1 / ( 2 t_ins/(Ch E_ins) + 2 JT/(Ch E_j) + SC_h/(Ch E_cbl) )
      E_r,turn = (2 JT E_j + 2 t_ins E_ins)/Cw + (SC_w/Cw) E_mid,r

- **Circonferenziale** (carico lungo la larghezza): le pareti superiore e
  inferiore sono in parallelo alla fascia centrale, nella quale isolante,
  pareti laterali del jacket e cavo sono in serie:

      E_mid,t = 1 / ( 2 t_ins/(Cw E_ins) + 2 JT/(Cw E_j) + SC_w/(Cw E_cbl) )
      E_t,turn = (2 JT E_j + 2 t_ins E_ins)/Ch + (SC_h/Ch) E_mid,t

- **Assiale**: media sulle aree, E_z,turn = (A_j E_j + A_cbl E_cbl + A_ins E_ins)/A_cella.

Nota: `cavo_stiffness` (in `size_cicc_cable`) usa una combinazione
diversa, che mescola le frazioni di larghezza e di altezza della cella
(per esempio 2 E_j JT/Cond_h nella rigidezza radiale, dove la parete
laterale occupa la frazione 2 JT/Cond_w); qui si usano le frazioni
geometriche corrette.

### 3.2 Anello

Con f_k la frazione della corda occupata dal componente k (pareti laterali
w_wall = corda − W − 2 GIT, isolanti laterali 2 GIT, WP W = n_turn Cond_w):

    radiale        E_r = Σ f_k E_r,k            (componenti in parallelo)
    circonferenz.  1/E_t = Σ f_k / E_t,k         (componenti in serie)
    assiale        E_z = Σ f_k E_z,k             (parallelo)

Gli anelli di acciaio (nose, piastra) sono isotropi con l'accoppiamento di
Poisson; gli anelli misti sono ortotropi senza accoppiamento di Poisson.
Rigidezze in deformazione piana generalizzata, ε = (ε_r, ε_θ, ε_z):

    σ = C ε,   isotropo: C = E/((1+ν)(1−2ν)) [1−ν ν ν; ν 1−ν ν; ν ν 1−ν]
               misto:    C = diag(E_r, E_t, E_z)

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

- La tensione circonferenziale è la stessa in tutti i componenti in serie
  dell'anello: il WP ha σ_θ dell'anello, e la deformazione circonferenziale
  del turn è ε_θ,turn = σ_θ / E_t,turn.
- Ogni parete del jacket ha una deformazione imposta nella direzione in cui
  è "in parallelo", una tensione imposta nella direzione in cui è "in serie",
  e ε_z = ε₀. Le altre due tensioni vengono dalla legge di Hooke 3D:

        ε_par = (σ_par − ν (σ_ser + σ_z))/E,   ε₀ = (σ_z − ν (σ_par + σ_ser))/E

  - pareti laterali: ε_par = ε_r (radiale), σ_ser = E_mid,t ε_θ,turn;
  - pareti superiore e inferiore: ε_par = ε_θ,turn, σ_ser = E_mid,r ε_r;
  - pareti laterali del case: ε_par = ε_r, σ_ser = σ_θ dell'anello.

- Pm del jacket del layer = massima Tresca tra parete laterale e parete
  superiore. Il nose e la piastra si linearizzano attraverso lo spessore
  come le linee di classificazione del FE (membrana = media, flessione =
  parte lineare di σ_r, σ_θ, σ_z).

L'effetto Poisson su σ_z è essenziale: con σ_z = E ε₀ il jacket risultava
~100 MPa troppo caricato assialmente (289 contro ~190 MPa sul design 10).

## 7. Varianti del percorso radiale del WP

- **`radial` = 'ring'** (default): il WP è incollato al case, il carico
  radiale si ripartisce con le pareti laterali (rigide) nell'anello.
- **`radial` = 'column'**: il WP è una colonna che scorre sulle pareti; la
  pressione radiale sul layer k è la forza di Lorentz dei layer sopra per
  unità di larghezza di ogni layer (dove un layer è più stretto, i turn in
  sbalzo appoggiano sulla spalla del case):

        p_k = Σ_{j<k} F_j/W_j + F_k/(2 W_k)

  con attrito sulle pareti (Janssen): la compressione circonferenziale σ_θ
  preme il WP sulle due pareti e ogni layer di altezza h_k cede al case fino a
  2 μ |σ_θ| h_k per unità di lunghezza (μ = `mu_case`, 0.2 come nel FE).
  Deformazione radiale del turn: ε_r = −p_k / E_r,turn.
- **`wall_radial`**: quota delle pareti laterali nella rigidezza radiale
  dell'anello (1 incollato, 0 scorrevole).

## 8. Risultati (02/10/2026)

### 8.1 Jacket, 23 run FE (Pm primario per layer)

| Variante | Rapporto medio modello/FE | rms | Massimo del design, min … max |
|---|---:|---:|---|
| anello, senza Poisson su σ_z | — | — | 15–20% sopra il FE |
| anello (default) | 0.86 | 18.3% | 0.62 … 1.08 |
| colonna senza attrito | 1.04 | 14.1% | 0.88 … 1.29 |
| colonna con attrito μ = 0.2 | 0.89 | 14.0% | 0.62 … 1.08 |
| surrogato (riferimento) | 1.0 | 7% | 0.87 … 1.15 |

**Risposta allo spessore del jacket** senza nessuna taratura: bench da +0.5
a +1.0 mm, FE −4.4% / −8.4%, modello −4.8% / −9.3%; d10 FE −5.5% /
−10.9%, modello −6.9% / −12.9%. Il meccanismo che conta per dimensionare
JT è quindi descritto bene.

### 8.2 Componenti delle pareti (design 7, turn centrali)

| Layer | FE parete lat. σ_r / σ_θ / σ_z | modello (anello) | FE parete sup. σ_θ / σ_r / σ_z | modello (anello) |
|---|---|---|---|---|
| L1 | +67 / −106 / 277 | +67 / −24 / 256 | −283 / 1 / 207 | −382 / 0 / 133 |
| L7 | −260 / −187 / 159 | −28 / −24 / 229 | −315 / −132 / 159 | −378 / −6 / 132 |
| L13 | −289 / −234 / 137 | −97 / −23 / 209 | −398 / −159 / 127 | −372 / −10 / 133 |

(FE = caso totale, con raffreddamento; MPa.) Il modello di cella in serie e
parallelo separa troppo i percorsi: nel FE le pareti laterali portano anche
compressione circonferenziale (fino a −234 MPa, modello −24) e le pareti
superiori anche compressione radiale (fino a −160, modello −10). Il
trasferimento passa per i raccordi arrotondati e per il cavo confinato.
L'anello (circonferenziale e assiale) è coerente con il FE; il punto debole
è la cella.

### 8.3 Case

Primo confronto, design 10: Pm del nose 526 contro 631 MPa del FE (−17%).
Confronto su tutti i design: in corso (run FE di riferimento).

## 9. Prossimo passo proposto

Sostituire il modello di cella in serie e parallelo con una **cella
elementare risolta numericamente** una volta per geometria (un turn con
raccordi, cavo, jacket e isolante, condizioni periodiche, tre deformazioni
medie unitarie): dà le rigidezze equivalenti esatte della cella e le
matrici di localizzazione (tensioni di membrana di ogni parete per unità di
deformazione media). È un calcolo piccolo (centinaia di elementi), uno per
grade, riutilizzabile; l'anello resta analitico.
