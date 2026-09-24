# Risposta alla revisione del codice (LEGGIMI_REVIEW, punti C00–C08)

Per ogni punto: esito, cosa è stato fatto nel codice, verifica eseguita.
I numeri sono del design 7 (EM_2D007) e del benchmark 2026, eseguiti in Octave 8.4.

| ID | Esito | In breve |
|---|---|---|
| C00 | chiarito | Il "surrogato" è un FEM 2D autonomo, non un modello ridotto: intestazione e documentazione corrette |
| C01 | **corretto** | Controlli di validità obbligatori; senza convergenza o con un controllo fallito il risultato è INVALID e la figura di merito non viene dichiarata soddisfatta |
| C02 | **implementato** | Attrito incrementale con scorrimento accumulato e sequenza di carico selezionabile (insieme / prima raffreddamento poi energizzazione) |
| C03 | **corretto** | Script di estrazione dei riferimenti FEM nel repository, soglie di accettazione esplicite, esito PASS/FAIL con errore |
| C04 | **implementato** | Due casi di carico: primario P (Lorentz + assiale) e totale P+Q (con raffreddamento); criteri Pm ≤ Sm, Pm+Pb ≤ 1.5 Sm su P, P+Q ≤ 3 Sm; Sm esplicito nell'Excel |
| C05 | **corretto** (un bug trovato) | Area con tolleranza, Jacobiano, residuo, equilibrio globale delle forze, equilibrio assiale, copertura SCL |
| C06 | aperto, decisione | Integrazione nello scan: resta da decidere (proposte sotto) |
| C07 | **corretto** | Il vault ora è confrontato con `S_amm_VT` |
| C08 | **corretto** | Sorgente arrotondata come nel FEM e ricerca del picco su tutto il contorno del cavo |

---

## C00 – Natura del modello

Corretto il testo: `wp_mech_surrogate` è un **modello a elementi finiti autonomo** (mesh, assemblaggio, contatti, soluzione), non un modello ridotto né una superficie di risposta. "Surrogato" indica solo che prende il posto della corsa ANSYS nella verifica di un punto di progetto. Il nome della funzione è stato mantenuto per non rompere le chiamate.

## C01 – Convergenza ed esito

- `solve_contacts` registra la convergenza **per ogni passo di carico** (`info.step`) e il numero di coppie che cambiano ancora stato. La tolleranza è ora un parametro esplicito (`surrogate_contact_tol`, default 0.5% delle coppie).
- Nuova struttura `out.checks` con `out.valid`. Se un controllo fallisce:
  - `fom.status` = `INVALID (…)`;
  - `fom.ok` = falso qualunque siano le utilizzazioni;
  - il plot e la console lo dichiarano;
  - `main_WP_TF_design` stampa un messaggio d'errore esplicito.

## C02 – Attrito e sequenza di carico

- L'attrito è ora **incrementale**: per ogni coppia si memorizza lo scorrimento accumulato g_s. Una coppia aderente ha T = −k_p (g_t − g_s) e scorre se |T| > μN; a fine passo g_s viene aggiornato (return mapping elastico-perfettamente-plastico).
- **Sequenza di carico** (`surrogate_load_sequence`):
  - 0 = tutti i carichi insieme, come le corse ANSYS usate per la validazione (default);
  - 1 = prima il raffreddamento, poi l'energizzazione in `surrogate_em_steps` passi.
- **Sensibilità sul design 7 (sequenza 1 contro sequenza 0):**
  - picco nel raccordo: differenze ≤ 3 MPa (1023 contro 1020);
  - Pm+Pb con tutti i carichi: ≤ 5 MPa (885 contro 882);
  - Pm: invariato.

  Per questa geometria la dipendenza dal percorso è quindi trascurabile. Per ogni caso nuovo resta disponibile l'opzione per verificarla.

## C03 – Riproducibilità dei benchmark

- `validation/tools/` contiene gli script Python con cui sono stati estratti i riferimenti FEM, con un README che spiega come ripeterli su una nuova corsa.
- `validate_mech_surrogate_2026.m` ha ora **bande di accettazione** esplicite:
  - εz entro ±2%;
  - picco per layer entro [0.85, 1.10];
  - Pm per layer entro [0.90, 1.10];
  - Pm+Pb per layer entro [0.85, 1.10];
  - SCL del case entro [0.90, 1.05] per nose e vault e [0.95, 1.05] per le pareti.

  Richiede anche che i controlli di validità passino. Termina con PASS o con un errore che elenca i fallimenti.
- **Esito della validazione con le bande: PASSATA** su entrambi i casi, con tutti i controlli di validità superati:

  | Rapporto modello FE / FEM | Benchmark | Design 7 |
  |---|---|---|
  | εz | +0.3% | 0.0% |
  | Picco per layer | 0.92–1.05 | 0.89–1.05 |
  | Pm per layer | 0.96–1.02 | 0.93–0.99 |
  | Pm+Pb per layer | 0.90–1.04 | 0.855–1.02 |
  | Nose / vault | 0.97 / 0.98 | 0.95 / 0.97 |
  | Piastra lato plasma | 1.14 | 1.06 |

  Il margine più stretto è Pm+Pb del layer 1 del design 7 (0.855 contro la soglia 0.85): il layer lato plasma, in flessione, è sottostimato del 15%.
- La forza assiale resta **prescritta dal FEM** per confrontare a parità di carico. Lo script stampa accanto il valore che lo strumento calcolerebbe dall'input.

## C04 – Classificazione delle tensioni

La linearizzazione da sola non classifica. Ora:

- si risolvono **due casi**:
  - **P**, primario: Lorentz + forza assiale, senza raffreddamento;
  - **P+Q**, totale: con il raffreddamento, che è il carico secondario;
- la figura di merito controlla:
  - Pm(P) ≤ Sm;
  - (Pm+Pb)(P) ≤ 1.5 Sm;
  - (Pm+Pb)(P+Q) ≤ 3 Sm;
  - su jacket e case;
- **Sm** è un parametro esplicito (`Sm_jacket`, `Sm_case`) e va **confermato** rispetto al codice di progetto e ai dati di materiale. Il default è pari a `S_amm_JT` / `S_amm_VT`, marcato "TO BE CONFIRMED" nell'Excel;
- con i contatti la separazione P / P+Q è approssimata: il caso primario è risolto con il proprio stato di contatto, senza il serraggio da raffreddamento, e questo è conservativo.

**Conseguenza sul design 7.** Con la classificazione, il primario dà Pm = 618 MPa (0.93 Sm), ma **Pm+Pb = 1107 MPa = 1.11 × 1.5 Sm**, nel layer 10 al turn di bordo.

- Senza raffreddamento il cavo non è serrato dal jacket e le pareti portano il carico di Lorentz in flessione.
- Il P+Q è 882 MPa (0.44 × 3 Sm).

Quindi **il design 7 non soddisfa il criterio primario di membrana + flessione**, se si accetta Sm = 667 MPa e si considera il raffreddamento come secondario. La valutazione precedente ("criteri soddisfatti") non lo classificava.

## C05 – Controlli numerici

| Controllo | Soglia | Design 7 |
|---|---|---|
| Area della mesh (isoparametrica) contro area del dominio | < 1·10⁻⁴ | 1.4·10⁻⁶ |
| Jacobiano in tutti i punti di Gauss | > 0 (altrimenti errore) | min 6·10⁻¹⁰ m² |
| Residuo della soluzione lineare | < 1·10⁻⁶ | 2·10⁻⁹ |
| Equilibrio globale (Lorentz + termico contro reazioni sui fianchi) | < 1·10⁻⁴ | 3.5·10⁻⁸ |
| Equilibrio assiale ∫σz dA contro T_bf | < 1·10⁻³ | 8·10⁻¹³ |
| Copertura dei punti sulle SCL del case | ≥ 0.98 | 1.00 |

Il controllo d'area ha trovato **un bug reale**: due triangoli spuri sovrapposti ai "ventagli" di filler, nel punto d'incontro tra due turn sul bordo del WP. Il baricentro cadeva esattamente sul confine tra due celle e il test con disuguaglianze strette non li escludeva. È stato corretto (intervalli chiusi). L'effetto sui risultati era trascurabile (circa 0.5 cm² di isolante in più su 2750 cm² di sezione), ma senza il controllo non sarebbe emerso.

## C06 – Integrazione nello scan

Resta aperto, ed è una scelta da fare.

- **Stato attuale:**
  - lo scan dimensiona ancora con il campo spalmato e con `SCF_transition_provisional` = 3.15;
  - la verifica FE è opzionale dopo la scelta della soluzione.
- **Proposte** (non ancora implementate, servono decisioni):
  1. dimensionare ogni grade sul picco di campo discreto, con un'iterazione tra dimensionamento e campo;
  2. usare il FE come filtro obbligatorio sulle soluzioni fattibili prima di salvarle come "accettate";
  3. tarare una stima analitica veloce sul FE per usarla nello scan al posto dello SCF provvisorio.

## C07 – Ammissibile del case

`search/scan_wp_designs.m`: la condizione di accettazione confronta ora `S_T_VT` con `S_amm_VT`; prima confrontava con `S_amm_JT`. Con i default attuali i due valori coincidono (667 MPa), quindi nessun risultato salvato finora cambia. Il bug si sarebbe manifestato con ammissibili diversi.

## C08 – Campo

- `wp_field_at_points` accetta il raggio di raccordo di ogni cavo e rappresenta la sezione arrotondata reale, a densità di corrente uniforme come nel FEM. La scomposizione è esatta in rettangoli: croce centrale più strisce a area conservata per ogni quarto di cerchio. La corrente totale è conservata: il rapporto del campo a grande distanza è 1.000001.
- `compute_discrete_field_profile` cerca il picco su **48 punti del contorno** arrotondato del cavo, dove cade il massimo dell'auto-campo, più la griglia interna 3×3.
- Risultato: 14.476 contro 14.482 T del FEM (design 7), 13.488 contro 13.489 T (benchmark).
- Anche le forze di Lorentz del modello meccanico usano ora la sorgente arrotondata.
