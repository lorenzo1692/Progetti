# Modello globale 3D a travi e shell (controparte MATLAB di STR_360)

Serve per la verifica a valle del dimensionamento 2D, in particolare per i
carichi fuori piano (forze toroidali da PF, torsione della gamba interna,
taglio nelle strutture di collegamento tra bobine). È la ricostruzione in
MATLAB del modello ANSYS `STR_360.dat` caricato il 29/09/2026.

## File

| File | Cosa fa |
|---|---|
| `coil3d/tf3d_global_model.m` | Costruisce il modello: tutte le n_TF bobine, travi lungo la linea media, shell tra bobine adiacenti (vault e OIS), supporti gravitazionali, carichi |
| `coil3d/tf3d_global_solve.m` | Soluzione statica lineare: travi di Timoshenko 3D, shell a 4 nodi (membrana con modi incompatibili QM6, piastra Mindlin MITC4, rigidezza di drilling), vincoli in assi locali |
| `validation/verify_tf3d_ansys_global.m` | Confronto con l'export ANSYS di `validation/ansys/POST_TF_BEAM_EXPORT.mac` |

## Modello (come STR_360)

- **Travi**: BEAM188 lungo la linea media di ogni bobina, sezione HREC
  0.675 × 0.675 m, parete 0.05 m (la cassa). Orientamento: asse locale z
  toroidale (il nodo K di ANSYS è spostato di 1 m in x per la bobina 1),
  quindi M_z è la flessione nel piano della bobina.
- **Shell**: SHELL181, t = 0.14 m, un elemento tra due bobine adiacenti per
  ogni tratto di linea media: gamba interna dritta (vault, nodi con R =
  R_min) e struttura di collegamento esterna (OIS, R tra 3.5 e 5 m, parte
  superiore e inferiore separate). Opzione `ois_wrap` = true riproduce
  l'elemento in più che STR_360 crea tra il primo e l'ultimo nodo di ogni
  gruppo OIS (non adiacenti lungo la bobina; `STR_360_improved` lo toglie).
- **Supporti**: per ogni bobina, spostamento toroidale e verticale nullo
  nel nodo più vicino a (R massimo, 0.8 · z minimo) tra i nodi con R ≤ 3 m.
- **Carichi**: forze nodali della bobina 1 (riga j di FL_TF_1 sul nodo j),
  ruotate su tutte le bobine; oppure un carico diverso per bobina.

## Verifica degli elementi (casi con soluzione analitica)

| Caso | Risultato |
|---|---|
| Mensola trave, carichi assiale, flessione nei due piani, torsione | esatta |
| Patch test di membrana su mesh distorta | spostamenti esatti, N uniforme |
| Flessione pura della membrana (una fila di elementi) | esatta (QM6) |
| Mensola di piastra, 8 elementi | −0.4% rispetto a Eulero-Bernoulli |

## Confronto con ANSYS STR_360 (carichi nel piano, run del 29/09/2026)

12 bobine × 162 nodi, 1944 travi, 924 shell, 11 664 gradi di libertà,
13 s in Octave.

| Grandezza | ANSYS | MATLAB |
|---|---|---|
| Forza assiale nella trave | 7.6 … 52.4 MN | 7.8 … 52.3 MN (scarto max 1%) |
| Momento nel piano | −2.18 … 3.00 MN·m | −2.17 … 2.92 MN·m |
| Taglio nel piano | −3.85 … 3.37 MN | −3.87 … 3.43 MN |
| Sforzo massimo in fibra della cassa (N/A + flessione) | 504 MPa | 500 MPa |
| Spostamento radiale / verticale | −5.02 … 0.39 / −0.39 … 10.59 mm | −4.93 … 0.32 / −0.18 … 10.41 mm |
| Reazione verticale al supporto, per bobina | 1.838 MN | 1.838 MN |
| Vault: N11 (circonferenziale) / N22 | −109 … −87 / 13.8 … 17.8 MN/m | entro 0.5% / 0.7% |
| Vault: M11, M22, Q23 | | entro 5–6% |
| Vault: intensità di Tresca max (superfici) | 883 MPa | 880 MPa |
| OIS: N22, M11, M22, Q23 | | entro 1–10% |
| OIS: N11 | fino a 6.9 MN/m | scarto max 1.1 MN/m (agli estremi della zona) |
| OIS: intensità di Tresca max | 107 / 133 MPa | 107 / 126 MPa |

Gli scarti massimi del momento e del taglio nel piano (12% e 7% del
massimo) sono locali, agli estremi della zona OIS (R ≈ 5 m), dove una sola
fila di shell attraverso la luce è sensibile alla formulazione
dell'elemento (SHELL181 integrazione ridotta contro QM6/MITC4).

Note sul confronto:
- l'export ANSYS dà come "centroide" degli shell il punto medio del lato
  sulla bobina 1, uguale per i due shell ai lati della bobina: il lato è
  riconosciuto dal numero di elemento (STR_360 crea prima la coppia 1-2,
  poi la coppia 12-1);
- STR_360 costruisce gli shell della coppia 12-1 e l'elemento di chiusura
  del vault a z = 0 con i nodi in ordine inverso: normale ribaltata, quindi
  M11, M22 e Q13 cambiano segno (M12 e Q23 no). Il confronto ne tiene conto.

## Carichi fuori piano

Con la componente toroidale di FL_TF_1 (risultante sulla bobina 1
−3.55 MN), MATLAB dà: torsione massima 3.65 MN·m sulla gamba interna
(R 1.12 m, z −1.8 m), flessione fuori piano 5.9 MN·m e taglio fuori piano
9.8 MN all'estremo della zona OIS (R ≈ 5 m), spostamento toroidale massimo
6.7 mm, N12 nel vault −18 … −9 MN/m. Manca ancora il confronto con ANSYS.

**Attenzione a `STR_360_improved.mac`**: applica la forza toroidale come FY
in coordinate nodali cilindriche (`F,N_CUR,FY,F_TOR(J)` con F_TOR = colonna
Fx globale di FL_TF_1). Sulla bobina 1 (θ = 90°) la direzione tangenziale
è −x globale, quindi la forza toroidale risulta applicata con il segno
opposto. Da solo il caso fuori piano cambia solo di segno; combinato con
altri carichi (gravità, guasti) il risultato è sbagliato. Correzione:
`F,N_CUR,FY,-F_TOR(J)` per la bobina 1, o in generale la proiezione della
forza globale ruotata sulla direzione tangenziale di ogni bobina.

## Collegamento alla run principale e dimensionamento accoppiato (03/10/2026)

### Il modello costruito dalla soluzione scelta

`coil3d/tf3d_global_from_design.m` costruisce il modello globale dalla
soluzione di MADE, invece che dai dati fissi di STR_360:

| Elemento | Da dove viene |
|---|---|
| Linea media | linea baricentrica della corrente del WP (`tf3d_design_geometry`), D bending free da r1 a r2, gamba interna dritta, `g3d_n_seg` = 160 tratti |
| Carichi | forze di Lorentz nodali dal campo 3D di tutte le bobine (`tf3d_centreline_loads`, come MAG_360), sezione del WP a 4 × 4 filamenti; PF/CS/plasma facoltativi (`opts.pf`) |
| Sezione della trave | sezione della gamba interna: case (settore tra il foro ad arco e la faccia lato plasma, meno la cava a inviluppo convesso come nel FE 2D) più i jacket spalmati sui layer (`tf3d_section_from_design`) |
| Shell del vault | spessore = acciaio che attraversa il piano medio della bobina: nose + piastra + pareti orizzontali dei jacket (2 JT per layer) |
| Shell OIS | `g3d_t_ois` = 0.14 m, come STR_360 |
| Zona OIS e supporto | alle stesse frazioni dell'escursione radiale della linea media di STR_360 (OIS tra R = 3.5 e 5 m, supporto tra i nodi con R ≤ 3 m, su 1.082–5.207 m) |

Il modello di STR_360 resta riprodotto come prima: la verifica
`verify_tf3d_ansys_global` dà gli stessi scarti (forza assiale 1%, Tresca
del vault 880 contro 883 MPa).

### Il collegamento con il dimensionamento 2D

Tutto il dimensionamento 2D carica la gamba interna con la forza verticale
della formula bending free:

    T_bf = 0.5 k_bf n_TF (N I)² μ0 / (2π),     k_bf = 0.5 ln(RTFo / RTFi)

cioè metà della forza verticale della metà superiore di una bobina. Il
modello globale dà la forza che la gamba interna porta davvero:

    F_inner(z) = N_trave(z) + N22_vault(z) · corda

Il loro rapporto è il **fattore di carico assiale**
k = max F_inner / T_bf (`axial_mode` = 1: con la flessione nel piano alla
fibra estrema). Il nuovo parametro `axial_load_factor` (default 1)
moltiplica T_bf in tutti i passi: scansione (tensione assiale del jacket,
`size_case_vault`), modello a strati e surrogato del case, surrogato del
jacket, `refine_jacket_fe`, FE 2D, export ANSYS (`AXIAL_LOAD_FACTOR`, letto
da `Solution_STR.f`).

### Dimensionamento accoppiato

`search/couple_global_sizing.m` chiude il ciclo:

1. modello globale della soluzione → k;
2. se k differisce di più dell'1% dal fattore con cui la soluzione è stata
   dimensionata, la stessa candidata viene dimensionata di nuovo con
   `axial_load_factor` = k (`search/resize_design_point.m`: stessa
   larghezza laterale del case, stessa combinazione iniziale di turn e
   layer, catena completa di cavo, jacket, nose, campo e controlli);
3. modello globale della soluzione ridimensionata (nuova sezione, nuovo
   baricentro della corrente), e di nuovo al punto 2 finché k non cambia.

La scansione salva ora la combinazione iniziale (`n_turns0`): il
ridimensionamento di una candidata richiede qualche secondo. Per i
risultati salvati prima, la candidata è cercata tra le combinazioni della
stessa larghezza laterale con lo stesso layout finale (circa 2 minuti).

### Risultati sulla macchina del template (scan completo del 03/10/2026)

| | #1 | #2 | #3 |
|---|---|---|---|
| Forza verticale della metà superiore, per bobina | 90.2 MN | 90.2 MN | 91.6 MN |
| T_bf della formula | 37.5 MN | | |
| Forza della gamba interna | 42.3 MN | | |
| k | **1.127** | 1.128 | 1.147 |
| ln(r2/r1) / ln(RTFo/RTFi) | 1.127 | 1.128 | 1.134 |
| Flessione nel piano della gamba interna (fibra estrema) | 22 MPa | 23 MPa | 19 MPa |
| Tresca degli shell del vault | 471 MPa | 479 MPa | 472 MPa |

- **Da dove viene k.** La gamba interna porta esattamente metà della forza
  della D a filamento sottile calcolata sui raggi della **linea
  baricentrica della corrente** (r1 = 1.092 m, r2 = 5.199 m). La formula
  di MADE usa invece i raggi della faccia lato plasma del case
  (RTFi = 1.259 m, RTFo = 5.032 m): ln(r2/r1)/ln(RTFo/RTFi) = 1.127. La
  forza verticale in più delle bobine discrete (+6.6% rispetto alla D
  continua) va alla gamba esterna.
- **Il carico non dipende dalla sezione.** k varia dello 0.3% tra la
  soluzione originale e quella ridimensionata: il ciclo converge in un
  solo ridimensionamento.
- **La flessione della gamba interna è piccola** (≈ 20 MPa alla fibra
  estrema, concentrata alla fine del tratto dritto): la D bending free fa
  il suo lavoro. Per questo `axial_mode` = 0 è il default.

**Dimensionamento accoppiato della #1** (k = 1.127, convergenza al primo
ridimensionamento):

| | dimensionata con T_bf | dimensionamento accoppiato |
|---|---|---|
| Ingombro radiale | 471.2 mm | 489.1 mm |
| Nose | 150 mm | 158 mm |
| Jacket per grade | 2.3 / 2.8 / 3.5 mm | 2.4 / 3.0 / 4.0 mm |
| Pm del case (surrogato) | 633 MPa | 634 MPa |
| Pm del jacket (surrogato) | 574 MPa | 586 MPa |

Con la forza verticale vera, la soluzione #1 dello scan cresce di 18 mm.
La #1 dimensionata con T_bf, verificata col FE 2D con k = 1, soddisfaceva
i criteri; con la forza vera la tensione assiale sale di circa 37 MPa
(da 292 a 329 MPa) e il case, al limite (636 MPa contro Sm = 667 MPa),
non li soddisferebbe più.

### Nella run principale

Passo 6b del main, dopo il passo 3D:

1. "Run the 3D beam + shell global model ...": modello globale della
   soluzione scelta, riassunto e grafici (`plot_tf3d_global`), salvataggio
   di `gm` nel .mat dei risultati;
2. se k differisce dal fattore usato: "Size this design again ... (coupled
   sizing)": ciclo accoppiato, grafici della convergenza, salvataggio di
   `cpl` (`cpl.row` è la soluzione accoppiata, `cpl.p` l'input con
   `axial_load_factor`);
3. verifica facoltativa col FE 2D della soluzione accoppiata (con la forza
   verticale vera);
4. l'export ANSYS (passo 7) scrive la soluzione accoppiata e
   `AXIAL_LOAD_FACTOR`.

Per scansionare tutta la macchina con la forza vera: impostare nel file di
input `axial_load_factor` al valore trovato (1.13 sulla macchina del
template). k è quasi costante tra le soluzioni (1.127–1.147).

### OIS: posizione e spessore scelti e controllati prima del run

Gli OIS si possono dare come pannelli `ois = [s_inizio s_fine t; ...]`
(frazioni della lunghezza della linea media: 0 = piano medio interno, ~0.5 =
piano medio esterno; `t` spessore in m). Senza `ois` si usa la finestra in R
`g3d_ois_frac` (parte alta e bassa, spessore `g3d_t_ois`). In entrambi i casi
`coil3d/tf3d_ois_zones.m` controlla **prima** del calcolo dei carichi:

- errore: fuori da [0 1], inizio ≥ fine, t ≤ 0, meno di 2 nodi, pannelli
  sovrapposti tra loro o al vault;
- avviso: casse che si toccano (larghezza libera `2R sin(π/n_TF) − W ≤ 0`),
  spessore maggiore della larghezza libera, PF più vicino di `ois_clearance`
  (0.1 m) se `pf` è dato.

Nella run principale (passo 6b) i pannelli vengono disegnati sulla linea media
con le tacche di s e si conferma con y; con n si scrivono nuovi pannelli
(`[0.30 0.40 0.10; 0.60 0.70 0.10]`), Invio = finestra di default, q = salta.
La scelta passa anche al dimensionamento accoppiato. Da riga di comando:

```matlab
gm = tf3d_global_from_design(row, p, struct('ois', [0.30 0.40 0.10; 0.60 0.70 0.10], ...
                                            'ois_plot', true, 'ois_confirm', true));
```

### Confronto fuori piano con ANSYS (serve una run APDL)

1. ANSYS: in `STR_360_improved.mac` correggere il segno della forza toroidale
   (bobina 1: `F,N_CUR,FY,-F_TOR(J)`, vedi sopra) e lanciare il SOLVE.
2. Nella stessa sessione: `/INPUT,'POST_TF_BEAM_EXPORT','mac'`
   (`validation/ansys/`); copiare `FL_TF_1.csv` nella stessa cartella.
3. MATLAB, da `MADE-2.0/matlab`:

   ```matlab
   addpath(genpath(pwd)); rmpath(genpath(fullfile(pwd,'legacy')));
   R0 = verify_tf3d_ansys_global('C:\run_ansys', 'verify_inplane', struct('with_toroidal', false, 'ois_wrap', false));
   R1 = verify_tf3d_ansys_global('C:\run_ansys', 'verify_oop',     struct('with_toroidal', true,  'ois_wrap', false));
   ```

   Risultati in `verify_oop/verify_tf3d_ansys_global.txt`. Attesi (MATLAB):
   torsione max 3.65 MN·m sulla gamba interna, flessione fuori piano 5.9 MN·m
   e taglio 9.8 MN a R ≈ 5 m, spostamento toroidale max 6.7 mm.

## Prossimi passi

1. Run ANSYS di `STR_360_improved` con la forza toroidale (segno corretto) e
   lo stesso export: confronto fuori piano.
2. Carichi PF nel modello costruito dalla soluzione (`opts.pf`), dai dati
   dei circuiti poloidali: torsione della gamba interna e taglio nell'OIS.
3. Stima di k già nella scansione, candidata per candidata, con i raggi
   della linea baricentrica (ln(r2/r1)/ln(RTFo/RTFi)), così lo scan parte
   già vicino alla soluzione accoppiata.
