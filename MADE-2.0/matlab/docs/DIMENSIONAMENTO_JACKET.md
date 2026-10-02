# Dimensionamento di jacket e case nello scan

Diario di lavoro: perché il dimensionamento va cambiato, cosa si fa, con
quali dati, e cosa resta aperto. I numeri di dettaglio stanno in
`validation/TF_FEM_benchmark_2026_findings.md`.

## 1. Obiettivo

Lo scan deve provare tutte le combinazioni (layout dei turn, larghezza
laterale del case, grading) e restituire con accuratezza la soluzione di
**minimo ingombro radiale** (Ri_ − Rk_) da cui partire. Vincolo: lo scan
resta snello, solo calcoli algebrici o surrogati. Il FE 2D
(`wp_mech_surrogate`, circa 10–20 minuti per design) si usa solo per
verificare la soluzione scelta, fuori dallo scan.

## 2. Come funziona oggi (versione del 29/09/2026)

1. Per ogni candidato, `size_cicc_cable` parte da `min_JT` e fa crescere
   JT in ogni grade finché la tensione di membrana sottile
   `p_rs·r_steel·dcr_jckt + S_z_JT` è sotto `S_amm_JT/safety_membrane`.
2. `size_case_vault` fa crescere il nose (DTF) finché la formula del vault
   (`S_z + S_c_VT`) arriva a `S_amm_VT`. Con `scf_model` = 0 nel ciclo entra
   anche la formula analitica del jacket con lo SCF provvisorio.
3. Con `scf_model` = 1 il jacket viene solo **verificato** dopo, con
   `jacket_stress_surrogate`: se non passa, il candidato è scartato.

Cosa ne risulta (scan del design 10, 31 design fattibili):

- **JT vale 2.1 mm in tutti i layer di tutti i design**, cioè `min_JT` più
  un passo: la formula di membrana sottile (413 MPa contro un limite di 513
  sul design 10) non governa mai. Il jacket in pratica non è dimensionato.
- Il vault è sempre al limite (S_T_VT 665–667 MPa), per costruzione.
- Il Pm+Pb del jacket dato dal surrogato sta tra 848 e 999 MPa, cioè appena
  sotto 1.5·Sm = 1000 MPa. Passano i layout che rientrano per caso; quelli
  che con un jacket un po' più spesso sarebbero buoni vengono scartati.
- Jacket e case non sono dimensionati insieme, anche se si influenzano: un
  jacket più spesso irrigidisce il WP e scarica il nose (vedi §3).

## 3. Sensitività al jacket (design 10, FE 2D, 01/10/2026)

Stesso layout, stesse aree di cavo, stesso nose (DTF 179 mm); sforzi
primari in MPa.

| JT per grade [mm] | Altezza WP | Jacket Pm | Jacket Pm+Pb | Surrogato Pm+Pb | Nose del case Pm |
|---|---:|---:|---:|---:|---:|
| 2.1 (scan) | 326 mm | 631 | 1075 | 962 | 631 |
| 2.6 | 351 mm | 596 | 995 | 846 | 613 |
| 3.1 | 377 mm | 562 | 925 | 759 | 597 |
| 2.1 / 2.6 / 3.1 | 350 mm | 581 | 948 | 833 | 611 |

- Il 50% di acciaio in più abbassa Pm+Pb del FE solo del 14%: il picco è la
  flessione locale del turn di bordo del layer 11 (14 turn) che sporge sul
  layer 12 (10 turn), poco sensibile a JT.
- Il surrogato scende del 21%: la sua dipendenza da JT è circa 1.5 volte
  troppo forte, quindi un margine fisso non è affidabile al variare di JT.
- Il jacket più spesso scarica anche il case: nose da 631 a 597 MPa a
  parità di spessore.

## 4. Piano

1. **Ritaratura del surrogato del jacket** (in corso, §5): aggiungere ai
   dati varianti degli stessi design a JT diversi, così che il fit separi
   l'effetto di JT da quello del layout; provare un termine per il gradino
   di larghezza; fissare il margine sull'errore residuo misurato.
2. **Jacket dimensionato per grade nello scan**: per ogni grade il JT
   minimo che rispetta Pm ≤ Sm e Pm+Pb ≤ 1.5·Sm (surrogato con margine),
   iterato con il nose (JT cambia altezza del WP, rigidezza, S_z). Solo
   formule: 2–3 passate per candidato.
3. **Grading come dimensione della ricerca**: numero di grade e campi di
   transizione (oggi fissi nell'input: `n_grades`, `grade_B_target_2/3`).
4. **Layout**: limitare o penalizzare il calo di turn tra layer adiacenti,
   che genera i picchi locali.
5. **Classifica** per minimo ingombro radiale; verifica FE 2D sulla prima o
   sulle prime 2–3 soluzioni, fuori dallo scan.

## 5. Passo 1: varianti a JT diverso per la taratura

Strumenti (`validation/`):

- `jacket_jt_variant.m`: lo stesso design con JT aumentato di ΔJT (scalare o
  per layer), con la stessa area di cavo e larghezza di cella del design di
  partenza e lo stesso spessore del nose; altezza del cavo e della cella
  ricalcolate dall'area, Rk_ spostato verso l'interno di quanto cresce il WP.
  Con ΔJT = 0 riproduce esattamente il design di partenza (verificato sui 6
  design base).
- `run_jacket_jt_calibration.m`: lancia i run FE (`jacket_surrogate_fe_data`)
  su più sessioni parallele, salva un file per run e salta quelli già fatti
  (si può riavviare dopo un'interruzione); con `worker = 0` aggiunge i
  risultati a `results/jacket_surrogate_calibration.csv`.
- `results/jacket_calibration_base_rows.mat`: i 6 design base.

Design base e varianti (13 run):

| Base | Layer | Turn per layer | JT base | Varianti |
|---|---:|---|---:|---|
| bench | 11 | 10 → 8 | 3.5 mm | +0.5, +1.0 mm |
| d7 | 13 | 8 | 3.0 mm | +0.5, +1.0 mm |
| d10 | 13 | 18 → 16 → 14 → 10 | 2.1 mm | +0.5, +1.0 mm, +0/+0.5/+1.0 per grade |
| s1_1 | 14 | 8 (nessun gradino) | da scan | +0.5, +1.0 mm |
| s1_13 | 25 | 12 → 10 → 8 | da scan | +0.5, +1.0 mm |
| s2_4 | 24 | 8 → 6 → 4 | da scan | +0.5, +1.0 mm |

Nel fit le varianti di uno stesso design (nome `<base>@jt+…`) vengono
escluse insieme nella verifica leave-one-design-out, altrimenti il test
vedrebbe il design che deve prevedere.

Esito dei run (01/10/2026): 12 su 13 completati.

- `s1_13@jt+1.0` non è calcolabile: con il jacket a +1 mm i layer profondi
  escono dal settore (gap toroidale −5.6 mm) e la mesh non si chiude (in
  Octave, che non ha la triangolazione vincolata, l'errore compare prima).
- `s1_1@jt+1.0` (gap 2.0 mm) e `s2_4@jt+1.0` (gap −0.6 mm) sono completi ma
  fuori dal settore: lo scan, che chiede `toroidal_gap` = 15 mm, avrebbe
  tolto dei turn. La parete laterale del case diventa un legamento sottile
  e lo sforzo FE **cresce** con JT (s2_4: 1165 → 1493 MPa). Restano nel CSV
  ma sono esclusi dal fit (soglia 5 mm). `jacket_jt_variant` ora restituisce
  il gap minimo (`min_toroidal_gap`).
- Correzione al FE trovata strada facendo: una linea di classificazione del
  case che non attraversa la mesh mandava in errore `trapz`; ora è segnata
  non valida dal controllo di copertura.

## 6. Nuova forma del surrogato

Dati: 392 layer, 23 run, 13 design base. Verifica leave-one-design-out (le
varianti di un design escluse insieme). Forme provate: esponente su JT
(0 … 1), termine per il gradino di larghezza, larghezza relativa del
layer, gradino del layer sopra, Cw/JT, Ch/Cw, numero di turn. Scelta:

    σ_k = S_z + (JT_k/3 mm)^0.2 · (a·σ_nom,k + b·σ_acc,k + c·σ_acc,k·s_k)
    s_k = (n_k − n_k+1)/n_k     (gradino di larghezza sotto il layer k)

| | a | b | c |
|---|---:|---:|---:|
| P<sub>m</sub> | 0.7532 | 0.2151 | 0.1983 |
| P<sub>m</sub>+P<sub>b</sub> | 1.3436 | 0.5489 | 0.5011 |

Le altre variabili non riducono l'errore in modo apprezzabile (rms per
layer da 10.2% a 9.4% nel caso migliore, margine richiesto quasi
invariato): lo scarto residuo dipende dal layout in un modo che le
grandezze algebriche dello scan non descrivono.

**Risposta a JT** (massimo Pm+Pb del design, rapporto con il design base):

| Design | JT | FE | Surrogato nuovo | Surrogato vecchio |
|---|---:|---:|---:|---:|
| bench +0.5 / +1.0 | 4.0 / 4.5 mm | 0.936 / 0.879 | 0.938 / 0.886 | 0.920 / 0.853 |
| d7 +0.5 / +1.0 | 3.5 / 4.0 mm | 0.917 / 0.850 | 0.924 / 0.861 | 0.900 / 0.820 |
| d10 +0.5 / +1.0 | 2.6 / 3.1 mm | 0.925 / 0.861 | 0.909 / 0.837 | 0.879 / 0.789 |
| s1_1 +0.5 | 4.0 mm | 0.934 | 0.936 | 0.916 |
| s1_13 +0.5 | 2.7 mm | 0.937 | 0.910 | 0.879 |
| s2_4 +0.5 | 4.2 mm | 0.908 | 0.941 | 0.921 |

Il surrogato nuovo segue l'effetto del jacket entro circa 3%; il vecchio
lo sovrastimava fino a 7 punti.

**Errore sul massimo del design** (FE/surrogato, 23 run):

| | mediana | 80° percentile | 90° percentile | max | rms per layer |
|---|---:|---:|---:|---:|---:|
| P<sub>m</sub> | 1.03 | 1.09 | 1.13 | 1.15 | 7.1% |
| P<sub>m</sub>+P<sub>b</sub> | 1.03 | 1.13 | 1.17 | 1.23 | 10.2% |

Per design base (P<sub>m</sub>+P<sub>b</sub>): sovrastimati quelli a layer
tutti uguali (d7 0.90, s1_1 0.91, s1_4 0.81, s1_7 0.80), sottostimati
quelli con gradini larghi (d10 1.19, s2_2 1.23, s2_4 1.12, s1_13 1.10).

## 7. Margine e dimensionamento nello scan

`jacket_margin` = 1.10 (default, nel file di input): copre 17 dei 23 run su
Pm+Pb e 20 su 23 su Pm. Un margine sul caso peggiore (1.23) farebbe
sovradimensionare quasi tutti i design di 10–20%, e la risposta debole di
Pm+Pb a JT (−14% con +50% di acciaio) lo trasformerebbe in molto acciaio e
ingombro radiale in più. La scelta è quindi: margine moderato nello scan
(classifica corretta, ingombro realistico) e verifica FE 2D della soluzione
scelta. Se la verifica non passa, la correzione si fa sul design stesso:
il rapporto FE/surrogato di quel design diventa il suo margine, si
ridimensiona JT con le formule e si conferma con un secondo run FE.

Nello scan (`scf_model` = 1, `search/scan_wp_designs.m`):

1. il candidato è dimensionato come prima (campo, cavo, JT minimo, nose);
2. `size_jacket_surrogate` trova per ogni layer il JT minimo con
   margin·Pm ≤ Sm e margin·(Pm+Pb) ≤ 1.5·Sm e assegna a ogni grade il
   massimo dei suoi layer (un conduttore per grade);
3. se il JT cresce, il candidato viene ridimensionato con quel JT come
   minimo (altezza del WP, gap toroidale, eventuale riduzione dei turn,
   nose, S_z) e il controllo si ripete, al massimo `jacket_max_passes`
   volte; poi la verifica del campo prosegue come prima;
4. se serve un JT oltre `JT_max` il candidato è scartato;
5. le soluzioni sono ordinate per ingombro radiale crescente (Ri_ − Rk_):
   la #1 è il punto di partenza.


## 8. Prova sulla macchina del design 10 (02/10/2026)

Scan ridotto con il nuovo dimensionamento (`scf_model` = 1,
`jacket_margin` = 1.10): layout con 18 turn nel primo layer, larghezza
laterale 0.0537 m come il design 10. 22 candidati, 6 fattibili, 6.5 minuti
in Octave. Classifica per ingombro radiale:

| # | Ingombro radiale | Layer | JT per grade [mm] | Surrogato Pm / Pm+Pb [MPa] |
|---|---:|---:|---|---|
| 1 | 542.4 mm | 13 (= design 10) | 2.1 / 2.1 / 2.4 | 558 / 908 |
| 2 | 548.4 mm | 14 | 2.1 / 2.4 / 2.5 | 550 / 891 |
| 3 | 592.8 mm | 18 | 2.1 / 2.3 / 2.4 | 551 / 896 |
| 4 | 604.4 mm | 19 | 2.1 / 2.2 | 556 / 898 |
| 5 | 637.8 mm | 20 | 2.1 / 2.2 | 558 / 899 |
| 6 | 692.6 mm | 25 | 2.1 | 556 / 899 |

La #1 ha lo stesso layout del design 10 (che con JT 2.1 mm ovunque non
passava il FE: Pm+Pb 1075 MPa); ora il grade 3 ha 2.4 mm e l'ingombro
cresce di 1.1 mm.

Verifica FE e correzione (`refine_jacket_fe`) sulla #1:

| Run FE | JT per grade [mm] | FE Pm | FE Pm+Pb | Surrogato Pm+Pb | Ingombro radiale |
|---|---|---:|---:|---:|---:|
| 1 | 2.1 / 2.1 / 2.4 | 626 | **1047** (L11) | 908 | 542.4 mm |
| 2 | 2.1 / 2.2 / 2.9 | 601 | **985** | 840 | 552.6 mm |

Al secondo run tutti i criteri sono soddisfatti (Pm ≤ 667, Pm+Pb ≤ 1000
MPa). Il design 10 corretto costa 10 mm di ingombro radiale rispetto alla
#1 dello scan. Questo layout è il caso peggiore della taratura per il
surrogato (FE/surrogato 1.15–1.19, gradino 14 → 10 turn): la correzione
FE è proprio il passo che serve per questi casi.

Nota sulla classifica: la #2 (548.4 mm, 14 layer) è sotto la #1 corretta
(552.6 mm) ma non è ancora verificata col FE; se anche lei richiedesse una
correzione l'ordine potrebbe cambiare. Per scegliere con certezza tra
soluzioni vicine (pochi mm) conviene verificare col FE le prime 2–3.
