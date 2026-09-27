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

Se il FEM non è valido, i plot restano (diagnostica) ma con banner rosso
"NOT VALIDATED" e messaggio che il FoM non va usato.

## 4. Verifica riproducibile

```matlab
p = read_machine_input('input/WP_TF_input_template.xlsx');
R = verify_ris_rect_fem(p, 'verify_out');   % report txt + PNG (sezione, GPS, zoom conduttore)
```

## 5. Limiti aperti

* ramo APDL RIS dell'export ANSYS non verificato in ANSYS;
* per RIS le sezioni di linearizzazione sono radiali: vicino agli angoli
  non sono normali alla faccia esterna;
* il caso RIS di verifica è un fixture sintetico dimensionato con
  `size_cicc_cable`, non l'uscita di una scansione;
* nel 3D il cavo RIS è rappresentato da 2×2 filamenti equivalenti (stesso
  momento del secondo ordine), non da un modello di autocampo del cavo
  tondo validato.
