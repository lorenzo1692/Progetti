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

## Come lanciarlo in MATLAB (senza ANSYS)

`coil3d/run_tf3d_global.m` usa il punto di progetto scelto in MADE:

```matlab
cd MADE-2.0/matlab; addpath(genpath(pwd)); rmpath(genpath(fullfile(pwd,'legacy')));
[row, p] = load_design_point('<tag>_results_<data>.xlsx', idx);   % la soluzione scelta
% 1) prima senza OIS: disegna la linea media con le tacche s/s_max per scegliere gli OIS
G = run_tf3d_global(row, p);
% 2) con gli OIS: [s_inizio s_fine spessore] per pannello, controllati e mostrati prima del run
G = run_tf3d_global(row, p, struct('ois', [0.30 0.40 0.10; 0.60 0.70 0.10]));
% con PF/CS (carichi fuori piano e controllo distanze): opts.pf = [rc zc dr dz I; ...]
```

- **Sezione delle travi** dal design (`coil3d/tf3d_beam_section.m`): cassa a
  cuneo tra `Rk_` e `Ri_` + WP strato per strato con il modulo longitudinale
  del turn (jacket, cavo, isolante), isolamento di massa; proprietà pesate sul
  modulo (E_ref = E_case); torsione con Bredt sulla cassa chiusa; aree di
  taglio dalle pareti. `section = 'hrec'` torna al quadrato cavo di STR_360.
- **Shell del vault**: spessore = naso `Rj_ - Rk_` (l'anello che regge la
  forza di centraggio nel dimensionamento 2D); `t_vault` per cambiarlo.
- **OIS** (`coil3d/tf3d_ois_zones.m`): posizione (frazione della lunghezza
  della linea media, 0 = piano medio interno, ~0.5 = piano medio esterno) e
  spessore di ogni pannello. Prima del run: errore se fuori da [0 1], sovrapposti
  tra loro o al vault, meno di 2 nodi; avviso se le casse si toccano (larghezza
  libera `2R sin(π/n_TF) − W ≤ 0`), se lo spessore supera la larghezza libera,
  se un PF è più vicino di `ois_clearance` (0.1 m). Poi grafico e conferma y/n
  (`confirm_ois = false` per saltarla).

## Come lanciare il confronto fuori piano

1. **ANSYS** – in `STR_360_improved.mac` correggere il segno della forza
   toroidale (vedi sopra): per la bobina 1 `F,N_CUR,FY,-F_TOR(J)`, oppure
   applicare la forza globale ruotata e proiettata sulla direzione
   tangenziale di ogni bobina. Lanciare il SOLVE.
2. Nella **stessa sessione** ANSYS, con la cartella di lavoro del run:
   `/INPUT,'POST_TF_BEAM_EXPORT','mac'` (file in `validation/ansys/`).
   Scrive `NODE_TF1_U.csv`, `BEAM_TF1_EL.csv`, `SHELL_TF1_EL.csv`,
   `REACT_TF.txt`, `SUM_TF.txt`. Copiare nella stessa cartella anche
   `FL_TF_1.csv` (i carichi usati da ANSYS).
3. **MATLAB**, dalla cartella `MADE-2.0/matlab` dello zip:

   ```matlab
   addpath(genpath(pwd)); rmpath(genpath(fullfile(pwd,'legacy')));
   % controllo: solo carichi nel piano, deve ridare la tabella sopra
   R0 = verify_tf3d_ansys_global('C:\percorso\run_ansys', 'verify_inplane', ...
        struct('with_toroidal', false, 'ois_wrap', false));
   % fuori piano: con la componente toroidale di FL_TF_1
   R1 = verify_tf3d_ansys_global('C:\percorso\run_ansys', 'verify_oop', ...
        struct('with_toroidal', true, 'ois_wrap', false));
   ```

   `ois_wrap = false` perché `STR_360_improved` non ha l'elemento OIS di
   chiusura (con lo `STR_360` originale: `true`). Se `n_TF`, sezione della
   cassa (`W`, `t_wall`) o spessore degli shell (`t_shell`) sono diversi dai
   default (12, 0.675 m, 0.05 m, 0.14 m) passarli nella stessa struct.
4. Risultati in `verify_oop/verify_tf3d_ansys_global.txt` + figure:
   per ogni grandezza (forze e momenti delle travi, N/M/Q degli shell,
   spostamenti, reazioni) range ANSYS, range MATLAB e scarto massimo.
   Attesi fuori piano (MATLAB): torsione max 3.65 MN·m sulla gamba interna,
   flessione fuori piano 5.9 MN·m e taglio 9.8 MN a R ≈ 5 m, spostamento
   toroidale max 6.7 mm.

## Prossimi passi

1. Run ANSYS di `STR_360_improved` con la forza toroidale (segno corretto) e
   lo stesso export: confronto fuori piano.
2. Collegamento al flusso di MADE: linea media dalla forma MATLAB, carichi
   dal modulo 3D (`tf3d_centreline_loads`), sezione della cassa dal
   dimensionamento 2D.
