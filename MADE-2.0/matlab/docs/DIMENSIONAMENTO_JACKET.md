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

Risultati: vedi §6 (da completare a fine run).
