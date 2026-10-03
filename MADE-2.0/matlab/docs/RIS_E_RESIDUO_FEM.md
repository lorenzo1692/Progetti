# RIS vs Rect e residuo del FEM surrogato

## 1. Perché RIS diventava Rect

`shape_cable` (200 = RIS, 201 = Rect) era usato **solo nel dimensionamento**
(`size_cicc_cable`: per RIS cella quadrata, JT = spessore minimo del jacket,
S_CICC con R_J = clamp(JT)). Non veniva salvato in `DATA`, quindi a valle
(plot sezione, campo discreto, FEM, export ANSYS, 3D) la geometria veniva
ricostruita sempre come rettangolo: cavo = Cond − 2JT − 2tins, che per una
cella RIS è un quadrato d×d invece di un cerchio di diametro d. L'export
ANSYS leggeva inoltre `p.shape_cable` corrente, non quello della soluzione.

## 2. Correzione

* `physics/wp_turn_geometry.m` (nuovo): unico punto in cui si ricava la
  geometria del turno da una soluzione. Forma presa da `row.shape_cable`,
  poi (solo se assente) da `p.shape_cable`; codici diversi da 200/201 →
  errore; RIS con cella non quadrata → errore (nessun fallback a Rect).
  * Rect: cavo SC_w×SC_h, raggio r = clamp(JT), jacket spessore costante JT.
  * RIS: cavo circolare d = a − 2JT = sqrt(4A/π) (a = Cond_w − 2 tins),
    jacket esterno quadrato raccordato R_J = clamp(JT), spessore variabile
    (JT a metà lato, massimo agli spigoli); isolante di turno offset del
    jacket (raggio R_J + tins); filler negli angoli della cella.
* `search/scan_wp_designs.m`: salva la colonna `shape_cable` in ogni riga.
* `postprocess/load_design_point.m`: se la colonna manca (file vecchi) usa
  il `p` salvato nel .mat della corsa; per xlsx la prende dall'input con
  warning `load_design_point:shape_from_input`. Un `p` passato in override
  non cambia la forma della soluzione (warning se differisce).
* Plot sezione / campo B (`plot_wp_section*.m`, `draw_turn_section.m`),
  campo discreto (`compute_discrete_field_profile.m`,
  `wp_field_at_points.m`), export ANSYS (`CABLE_TYPE 'RIS'`, ramo APDL
  non verificato → warning `export_ansys_input:ris_unverified`), modulo 3D
  (`tf3d_design_geometry.m`, `tf3d_coil_filaments.m`: per RIS i punti di
  Gauss usano il quadrato equivalente di lato d·√3/2, stesso momento del
  secondo ordine del cerchio d²/16).
* FEM `physics/wp_mech_surrogate.m`:
  * mesh RIS: stazioni sul contorno esterno del jacket proiettate
    radialmente sul cerchio del cavo; strati del jacket per interpolazione
    lineare fra cerchio e quadrato raccordato (spessore variabile); isolante
    come offset; 2 elementi per raccordo esterno; punti interni per la
    triangolazione del cavo (T6);
  * contatti cavo/jacket sui nodi del cerchio (nodi duplicati, Coulomb);
    il lato del cavo verso cui spinge la forza di Lorentz del turn
    (± `surrogate_ris_bond_angle`, default 45° = un lato del quadrato) è
    sempre chiuso (bilaterale in normale) e la sua coppia centrale è legata
    anche in tangenziale: un cerchio in un foro circolare non ha rigidezza
    a rotazione e, se tutte le coppie si aprono (raffreddamento prima
    dell'energizzazione), nemmeno a traslazione. La trazione trasmessa dal
    lato legato è riportata (`bond_tension_fraction`) e deve restare
    piccola perché l'ipotesi "lato premuto" valga; `surrogate_ris_bond = 0`
    torna al contatto unilaterale ovunque;
  * area del cavo per la densità di corrente = area del cerchio, verificata
    sulla mesh (`cable_area_ok`, errore relativo < 1e-3);
  * carichi EM e risultanti (`out.load_resultant`), sezioni di
    linearizzazione e post-processing sul contorno del jacket reale.

## 3. Residuo (FAILED: residual 4.1e-6)

Non è stata alzata la tolleranza (1e-6) né toccato il criterio del FoM.
Diagnosi: il residuo della prima soluzione Cholesky non scalata era
dominato dal condizionamento (penalty 1e3·max(diag K) su una K già
penalizzata dai tie → penalty dei contatti amplificata, diagonale con
scale molto diverse, incluso eps_z). Correzioni:

1. penalty di riferimento k_ref = max diagonale di spostamento della K
   elastica (prima dei tie, escluso eps_z); tie e contatti = 1e3·k_ref;
2. Cholesky sparso su matrice scalata di Jacobi (simmetrica) + fino a 2
   passi di raffinamento iterativo; il residuo riportato è quello di U
   sul sistema **non scalato**: ‖K U − F‖/‖F‖;
3. residuo e bilancio controllati su caso totale **e** primario (max);
4. nuovo controllo di consistenza dello stato di contatto sulla U finale
   (`contact_state_ok`): forza delle coppie che cambierebbero stato
   (normale e attrito) / forza normale totale < 2e-2. Il bilancio globale
   da solo non prova l'equilibrio locale: per questo si controllano
   residuo, bilancio e stato di contatto separatamente.

5. la convergenza dell'iterazione di contatto richiede anche la coerenza
   di stato sotto `contact_violation_tol` (2e-2, la stessa soglia del
   controllo): nel benchmark A il caso primario si fermava con coppie in
   bilico che portavano il 5% della forza normale (residuo 1.5e-9, ma
   stato non coerente); ora continua a iterare (entro `contact_maxit`) e
   arriva a 2.7e-5. Validazione ANSYS A + design 7: PASSED, rapporti
   invariati.

### Penalty o Lagrangiano aumentato

Penalty: N = -kp*gn, penetrazione gn = N/kp. Lagrangiano aumentato:
N = λ - kp*gn con λ aggiornato (Uzawa) fino a penetrazione nulla entro
tolleranza, vincolo esatto anche con kp basso (matrice meglio
condizionata), ma più soluzioni per iterazione di stato e moltiplicatori
anche tangenziali per l'attrito. Qui la penalty è sufficiente perché il
difetto che l'AL curerebbe (condizionamento) è risolto dallo scaling di
Jacobi (residuo ~1e-9 con kp = 1e3 k_ref) e l'errore della penalty è
misurato: `penetration_rel` = penetrazione massima / JT minimo, riportata
in console e nel report di verifica. L'AL conviene se la penetrazione
supera ~1% di JT o se serve kp più basso; si può aggiungere come opzione
aggiornando λ a ogni iterazione di stato.

Se il FEM non è valido, i plot restano (diagnostica) ma con banner rosso
"NOT VALIDATED" e messaggio che il FoM non va usato.

## 4. Campo di picco nella scansione

Il campo con cui la scansione dimensiona i grade era Ampère × `corr_B_WP`
(costante, 1.05 nel template) con profilo lineare. Contro il modello
discreto validato su ANSYS (0.1%) sottostima il picco di ~1 T già per un
WP largo (design 7: 13.48 contro 14.48 T) e di 1.7-2 T per D/W > 1.5; il
rapporto picco/Ampère va da 1.12 a 1.20, dipende da D/W, da Iop
(auto-campo) e dalla gradazione, e la geometria è nota solo dopo il
dimensionamento: nessun fattore costante o funzione del solo rapporto
d'aspetto lo corregge entro pochi decimi di tesla.

`physics/wp_peak_field_fast.m` calcola il picco per layer con lo stesso
modello (turn di bordo, vicino al bordo e centrali, 48 punti sul
contorno; rettangolo equivalente lontano, scomposizione esatta entro 2.5
celle), entro 0.09 T dal profilo completo, ~0.4 s in Octave.
`scan_wp_designs` (p.field_model = 'discrete', default) lo chiama sui
candidati che passano tutti i controlli e ridimensiona i grade al picco
reale finché |picco - campo di dimensionamento| <= p.field_tol (0.05 T),
al massimo p.field_max_iter (6) passate; nuove colonne B_peak,
B_peak_layers, field_iter. p.field_model = 'smeared' riproduce il
comportamento precedente.

Esito sul design 7 (scansione mirata, stessa macchina): con il campo
reale il grade 1 va dimensionato a 14.48 T e si estende ai layer 1-7; il
cavo Nb3Sn diventa alto 48.5 mm in una cella larga 43 mm (rapporto 0.89 <
0.99): il candidato non è più fattibile. Era accettato solo perché il
campo era sottostimato.

Scansione completa della macchina del template, tre larghezze di case
con soluzioni (0.134/0.154/0.174 m), dopo la correzione di `cicc` sotto:
campo spalmato 100 soluzioni (26/47/27), picco discreto 42 (10/25/7);
campo di dimensionamento del grade 1 13.76-14.20 T invece di 13.48 T;
area del cavo del grade 1 +15...+35%; 3-6 passate, tutte convergenti.

Correzione in `cicc.m`: la ricerca del rame si fermava sul punto con THS
più vicino al limite in valore assoluto, fino a 5 K sopra; il controllo
THS della scansione scartava il candidato per un arrotondamento (14 dei
20 progetti ammessi dal campo spalmato, con THS 251-255 K). Ora tiene il
più vicino al limite da sotto (più rame, mai meno). Col campo spalmato le
soluzioni sulle stesse tre larghezze passano da 20 a 100.

### Pressione magnetica dal campo del modello scelto

`p_rs` (dimensionamento del jacket, tensione radiale, vault) ora viene dal
campo con cui è dimensionato il grade 1: picco tarato o discreto con
`field_model` 1/2, Ampère × `corr_B_WP` con `field_model` 0 (invariato).
Con il campo tarato cresce di (k/1.05)^2, fino a circa +11% sui WP stretti.
La formula di Ampère (`2*pi*RTFi - dr_plasma_side`) resta com'è.

### Rapporto d'aspetto massimo della cella

`max_cable_aspect_ratio` (WP dimensioning, default 2): la scansione scarta
le celle con Cond_w/Cond_h sopra questo valore.

## 4b. Surrogato veloce della tensione del jacket (`scf_model` = 1)

Il metodo attuale (`scf_model` = 0, default) stima la tensione del jacket
con σ_nom × SCF di tabella (su Cond_w/SC_w) × `SCF_transition_provisional`
(3.15, sui layer vicini a un cambio di grade) × `dcr_WP_rad`. Contro il FE 2D
validato su ANSYS sottostima il Pm+Pb primario massimo del design fino al
62%: non vede la crescita della tensione con la profondità (circa +70% dal
layer lato plasma a quelli profondi) e `dcr_WP_rad` introduce una
variazione di 2 volte tra design con tensioni FE simili.

`physics/jacket_stress_surrogate.m` stima per ogni layer i Pm e Pm+Pb
primari del jacket come nel FE:

    stress_k = S_z + a*sigma_nom_k + b*sigma_acc_k
    Pm:    a = 0.711, b = 0.244      Pm+Pb: a = 1.285, b = 0.601

con sigma_nom = p_rs*r_steel*dcr_jckt (la formula attuale senza SCF) e
sigma_acc dalla pressione radiale accumulata vera (somma di n*Iop*B dei
layer sopra, per unità di larghezza del WP). La scansione li verifica come
la figura di merito del FE (Pm ≤ Sm_jacket, Pm+Pb ≤ 1.5 Sm_jacket) e, in
questa modalità, il nose è dimensionato solo dai criteri del vault. Nuove
colonne `JT_Pm`, `JT_PmPb` [MPa].

Taratura: 242 layer di 13 design a cavo rettangolare calcolati con il FE
(design 7, benchmark, design 10, 10 soluzioni della scansione tarata, W 304–344 mm,
14–29 layer, 21–66 kA). Verifica leave-one-design-out, errore sul massimo del
design: Pm+Pb −25…+13% (rms per layer 12%), Pm −16…+12% (rms 7%); formula
attuale: Pm+Pb −62…+6%. Il termine di transizione tra grade non è
significativo (il FE non mostra il salto che il 3.15 assume). L'errore
peggiore è sui WP stretti con gli ultimi layer ristretti (picchi locali di
flessione ai gradini di larghezza). Non tarato per RIS: con `scf_model` = 1
e `shape_cable` = 200 la scansione si ferma con un errore. È un modello
preliminare di screening: la soluzione scelta va verificata con il FE.

Riproducibilità: `validation/jacket_surrogate_fe_data.m` (dati FE per
layer), `validation/results/jacket_surrogate_calibration.csv`,
`validation/tools/fit_jacket_surrogate.py` (coefficienti ed errori).

Nota importante dai dati: il Pm+Pb primario del jacket nel FE va da 991 a
1369 MPa; 11 design su 12 superano 1.5 Sm = 1000 MPa, come già il design 7.

## 5. Verifica riproducibile

```matlab
p = read_machine_input('input/WP_TF_input_template.xlsx');
R = verify_ris_rect_fem(p, 'verify_out');   % report txt + PNG (sezione, GPS, zoom conduttore)
```

## 6. Limiti aperti

* ramo APDL RIS dell'export ANSYS non verificato in ANSYS;
* per RIS le sezioni di linearizzazione sono radiali: vicino agli angoli
  non sono normali alla faccia esterna;
* il caso RIS di verifica è un fixture sintetico dimensionato con
  `size_cicc_cable`, non l'uscita di una scansione;
* nel 3D il cavo RIS è rappresentato da 2×2 filamenti equivalenti (stesso
  momento del secondo ordine), non da un modello di autocampo del cavo
  tondo validato.
