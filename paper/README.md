# Paper MADE 2.0 — Nuclear Fusion

Bozza del paper *"MADE 2.0: a Magnet Design Explorer to derive superconducting
magnet systems from a target magnetic configuration for tokamaks, mirrors and
stellarators"* per **Nuclear Fusion** (IOP/IAEA).

## Messaggio centrale

1. **La configurazione magnetica è l'input, i magneti sono l'output.** Il fisico
   fissa B0, R0, A, ripple, flux swing, correnti PF, profilo B(z) del mirror;
   MADE restituisce conduttore, winding pack e struttura che la realizzano.
2. **Esplorare, non ottimizzare.** Si tengono tutti i punti di progetto
   fattibili → mappe di trade-off, la scelta resta al progettista.
3. **Stessi moduli, macchine diverse.** Conduttore + WP + struttura sono
   comuni a TF, CS, PF; il mirror riusa CS/PF (solenoidi + anelli);
   lo stellarator riusa il conduttore del TF *senza wedge* + modello globale.

## Struttura (`main.tex` + `sections/`)

| # | Sezione | Contenuto | Stato |
|---|---------|-----------|-------|
| 0 | Abstract | bozza completa, manca il numero della verifica | bozza |
| 1 | Introduction | piano per paragrafi | da scrivere |
| 2 | Design philosophy | config → magneti, explore vs optimise, architettura modulare (Fig. 1–2, Tab. 1) | parziale |
| 3 | Conductor module | selezione SC (NbTi/Nb3Sn/REBCO/ibrido), hot spot, giacca, rigidezze | bozza con equazioni |
| 4 | TF coils | R1/R2 da ripple, Princeton-D, grading, case/vault, scansione | bozza con equazioni |
| 5 | CS e PF | **codice non presente nel repo** | da scrivere |
| 6 | Mirror | riuso moduli CS/PF, flusso di progetto | schema |
| 7 | Stellarator + modello globale beam–shell | livello locale (TF senza wedge) + livello globale; primo caso: TF a D | schema dettagliato |
| 8 | Verification | confronto con design di riferimento | da fare |
| 9 | Applications | mappe di design space (DEMO/VNS/CEFTR…) | da fare |
| 10 | Conclusions | — | ultima |
| A | Appendice | rigidezze equivalenti del conduttore | bozza |

Tutti i punti aperti sono marcati `\todo{...}` (rosso nel PDF).

## Ruolo del modello globale beam + shell

È il pezzo che manca e che sblocca due cose insieme:

- **TF a D**: elimina l'ipotesi Princeton-D (solo membrana), dà flessione e
  carichi fuori piano con i PF → serve anche come **caso di verifica** contro FE
  di riferimento prima degli stellarator.
- **Stellarator**: bobine non planari senza vault a cuneo; il carico va in
  case + struttura intercoil. Beam curvi 3D lungo la centerline (proprietà di
  sezione omogeneizzate dal WP — `Ke_rad`, `Ke_tor` già calcolati in
  `size_cicc_cable.m`), shell per la struttura intercoil, carichi da
  Biot–Savart. Le forze di sezione tornano al livello locale per le verifiche.

Scelte da fare: frame locale (non Frenet, usare il frame finite-build della
bobina), solutore FE proprio in MATLAB vs codice FE esterno pilotato da MADE.

## Note emerse leggendo il codice (da sistemare prima dei numeri del paper)

- `compute_operating_params.m`: picco di campo con denominatore
  `2*pi*RTFi - dr_plasma_side`; probabilmente dovrebbe essere
  `2*pi*(RTFi - dr_plasma_side)`. Da verificare.
- `cicc.m`: commento "NbTi sotto 5 T" ma la soglia nel codice è `B < 6`.
- Margine operativo (I/Ic o temperatura) implicito: va dichiarato nel paper.
- `Tau_discharge = max(..., 4)`: giustificare il minimo di 4 s e
  l'espressione di `Tau_discharge1` (vincolo vacuum vessel).

## Domande aperte

1. Dove sono i solutori **CS e PF**? Vanno importati nel repo per scrivere la
   sez. 5 coerente con il codice.
2. Esiste una **MADE 1.0** pubblicata da citare?
3. Casi di **verifica** e macchine pubblicabili (DEMO, DTT, VNS, CEFTR?).
4. Caso mirror di riferimento (es. mirror HTS ad alto campo).
5. Lista autori/affiliazioni e funding (EUROfusion?).

## Compilazione

Serve `iopart.cls` dal template autori IOP (Nuclear Fusion). In alternativa
usare la riga di fallback `article` in `main.tex`.

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```
