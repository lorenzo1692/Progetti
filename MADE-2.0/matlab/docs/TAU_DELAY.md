# Criterio per il tau delay (2026-10-03)

## Come si sceglie (file di input)

Nel foglio `Input`, categoria *Thermal allowables*:

| Nome | Valore | Effetto |
|---|---|---|
| `Tau_delay_mode` | 0 | ritardo fisso `Tau_delay` per tutti i grade |
| `Tau_delay_mode` | 1 | criterio ITER: dump 2 s dopo che la zona normale ha raggiunto 0.1 V (`t_det` dalla fisica del cavo, sotto); i grade HTS usano `Tau_delay` (ITER non ha dati HTS) |
| `Tau_delay` | [s] | ritardo fisso (modo 0, e HTS nel modo 1) |
| `v_quench_LTS` | [m/s] | velocità di propagazione della zona normale nel CICC LTS (solo modo 1) |

Fonte dei valori ITER TF (0.1 V + 2 s): R. Zanino et al., *Quench analysis
of an ITER TF coil* (analisi di quench della bobina TF ITER: propagazione
fino a t(0.1 V) + 2 s). Altre fonti ITER: soglia 0.2 V su un pancake per il
TF; CS: soglia ±0.55 V e holding time 1.5 s (M. Coatanéa et al., IEEE TAS
25(3) 2015, DOI 10.1109/TASC.2015.2390296). **Da verificare sul testo
completo / ITER DDD prima della pubblicazione.**

Esempio (12 T, 60 kA, τ_dis 20 s, v_q 5 m/s): modo 1 → τ_delay 2.16 s,
516 strand di Cu; modo 0 con 1 s → 477; con 3 s → 545.

## Studio originale (proposta)

## Cosa rappresenta

In `heat_balance_cicc_ode` il punto caldo parte da Tc a t = 0 (innesco del
quench) e resta a corrente piena fino a `Tau_delay`, poi la corrente scende
con `Tau_discharge`. `Tau_delay` è quindi il tempo tra l'innesco locale e
l'inizio della scarica. Oggi vale 1 s fisso (`cicc_params.m`); in FED 2023 si
usavano 3 s, nel CS 0.5 s.

## Criterio proposto

`Tau_delay = t_det + t_hold + t_act`

| Termine | Significato | Da dove viene |
|---|---|---|
| `t_det` | tempo perché la zona normale sviluppi la tensione di soglia `V_th` | fisica del cavo (sotto) |
| `t_hold` | tempo di convalida del segnale (filtra disturbi induttivi: disruzioni, rampe PF) | sistema di rivelazione |
| `t_act` | apertura dell'interruttore / inserzione della resistenza di scarica | circuito di protezione |

Rivelazione in tensione: la zona normale cresce con due fronti a velocità
`v_q`, lunghezza `2 v_q t`, tensione `V(t) = 2 v_q t · I · ρ_Cu(T_ref, B) / A_Cu`.
Ponendo `V = V_th`:

    t_det = V_th · A_Cu / (2 · v_q · I · ρ_Cu(T_ref, B))

- `ρ_Cu` con lo stesso fit (RRR, magnetoresistenza) dell'ODE, alla temperatura
  iniziale (conservativo: la zona normale è più calda, quindi più resistiva).
- **Accoppiamento con il rame**: più rame → meno tensione → rivelazione più
  lenta. Il criterio è calcolato dentro la ricerca del rame di `cicc.m`, quindi
  il rame trovato è coerente con il ritardo che lo protegge.

Implementazione: `conductor/quench_delay_time.m`; attivazione con
`cp.Tau_delay_mode = 'detection'` in `conductor/cicc_params.m` (default
ancora `'fixed'`, nessun risultato cambia finché non si decide).
Studio: `validation/study_tau_delay.m`.

## Risultati (Iop 60 kA, τ_dis 20 s; LTS Nb3Sn a 12 T, HTS REBCO a 18 T)

Valori provvisori: `V_th` 0.1 V, `t_hold` 1 s, `t_act` 0.5 s.

| Caso | N_Cu | t_det [s] | Tau_delay [s] | Cu [mm²] |
|---|---|---|---|---|
| LTS fisso 1 s | 477 | – | 1.00 | 533 |
| LTS fisso 3 s | 545 | – | 3.00 | 586 |
| LTS rivelazione, v_q 1 m/s | 521 | 0.81 | 2.31 | 567 |
| LTS rivelazione, v_q 5 m/s | 500 | 0.16 | 1.66 | 551 |
| LTS rivelazione, v_q 20 m/s | 496 | 0.04 | 1.54 | 548 |
| HTS fisso 1 s | 837 | – | 1.00 | 657 |
| HTS rivelazione, v_q 0.2 m/s | 1014 | 3.95 | 5.45 | 796 |
| HTS rivelazione, v_q 0.05 m/s | 1632 | 25.5 | 27.0 | 1282 |
| HTS rivelazione, v_q 0.01 m/s | 6134 | 478 | 480 | 4818 |

## Cosa ne esce

1. **LTS**: `t_det` è piccolo (0.04–0.8 s) rispetto a `t_hold + t_act`; il
   ritardo è dominato dal sistema di protezione. L'effetto sul rame è ±5%.
   Il criterio è robusto: l'incertezza su `v_q` pesa poco.
2. **HTS**: con rivelazione in tensione `t_det` domina e diverge (più rame →
   rivelazione più lenta → più rame). Sotto ~0.1 m/s il cavo non è
   proteggibile in tensione. Per l'HTS il ritardo va quindi fissato dalla
   tecnologia di rivelazione (fibre ottiche, co-wound, ecc.):
   `cp.quench.t_det_HTS`.
3. Il paper può presentarlo così: per LTS `Tau_delay` derivato dalla fisica del
   cavo + requisiti del sistema di protezione; per HTS MADE restituisce il
   **massimo ritardo ammissibile** come requisito per la rivelazione.

## Da decidere / fonti da fissare

- `V_th`, `t_hold`, `t_act`: valori di progetto (ITER DDD magneti / criteri
  DEMO). Quelli nel codice sono segnaposto.
- `v_q` per CICC LTS: la propagazione è dominata dall'elio (espansione) e
  accelera nel tempo; un valore costante è una semplificazione. Possibile
  passo successivo: formula di Dresner/Shajii–Freidberg al posto di `v_q`
  costante.
- Possibile output aggiuntivo per ogni punto di progetto: il ritardo massimo
  ammissibile con il rame scelto (inverso del problema), utile soprattutto
  per l'HTS.
