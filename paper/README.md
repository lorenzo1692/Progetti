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
| 5 | CS e PF | da FED 2023, TAS 2024 e codice CS/PFC | bozza |
| 6 | Mirror | riuso moduli CS/PF, flusso di progetto | schema |
| 7 | Stellarator + modello globale beam–shell | modello già nel codice, verificato vs ANSYS nel piano (Tab.) | bozza |
| 8 | Verification | DTT, ITER, CEFTR | tabella input impostata |
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

Il solutore esiste già in MATLAB (`coil3d/tf3d_global_model.m` +
`tf3d_global_solve.m`, controparte di STR_360) ed è verificato contro ANSYS
nel piano. Mancano: confronto fuori piano, sezioni collegate al
dimensionamento 2D, frame locale per bobine non planari (finite-build, non
Frenet), carichi bobina per bobina senza simmetria di rotazione.

## Codice di riferimento (branch)

| Parte | Branch | Cartella |
|---|---|---|
| TF (versione evoluta, con correzioni), conduttore, modulo 3D, modello globale beam+shell | `ccr-165c7b01-g5w4kx` (= `claude/code-improvement-3cmng0` + fix) | `MADE-2.0/matlab/` (`physics/`, `conductor/`, `coil3d/`, `validation/`) |
| CS e PF | `claude/tf-cs-pfc-folders-ih6f07` | `MADE-2.0/matlab/CS`, `MADE-2.0/matlab/PFC` |
| Vecchi script per macchina (CEFTR 2025, VNS, DEMO…) | stesso branch | `MADE-2.0/matlab/legacy/scripts` |

## Stato dei problemi nel codice (aggiornato 2026-10-03)

Il branch `ccr-165c7b01-g5w4kx` contiene ora il TF evoluto
(`claude/code-improvement-3cmng0` mergiato) con le correzioni:

| # | Problema | Stato |
|---|---|---|
| 1 | Campo di Ampère `2*pi*RTFi - dr` | **Corretto**: `2*pi*(RTFi - dr)`; `B_PHI_0` dal campo a `RTFi` (ora ≈ B0: 5.706 vs 5.7 T sul template). Corretti anche `legacy/scripts` CEFTR/VNS. Template VNS: B_Amp +1.4% |
| 2 | Soglie conduttore | TF risolto (parametro); CS/PFC hanno ancora 6.2/14.5 T, d 0.82 mm, VF 0.8 → da unificare |
| 3 | Minimo 4 s su τ scarica | **Tolto** (refuso): `τ = max(τ_VV, L·Iop/V_max)` caso per caso (template: τ_VV = 1.63 s) |
| 4 | τ delay | **Proposta di criterio** in `MADE-2.0/matlab/docs/TAU_DELAY.md` (`conductor/quench_delay_time.m`, opzionale) |
| 5 | Margine operativo implicito | da dichiarare |

## Lavoro MADE già pubblicato (base del paper 2.0)

- MADE 1.0 — FED 193 (2023) 113659: macro-blocchi I (TF) e II (CS), caso DEMO.
- Macro-block III — IEEE TAS 2024, DOI 10.1109/TASC.2024.3362753: forma TF ottimizzata + PF, caso ITER.
- CS di DTT progettato con MADE (IEEE TAS 2020, 2021).

Novità 2.0 da rivendicare: codice modulare per famiglia di bobine; campo TF
discreto calibrato (validato ANSYS); surrogati meccanici FE; FEM
assialsimmetrico CS/PF interno; modulo 3D + modello globale beam+shell
verificato contro ANSYS (STR_360); mirror; percorso verso stellarator.

## Macchine per la verifica

DTT, ITER, CEFTR. ITER ha già il caso PF in Macro-block III (Tab. II) e il
DDD pubblico; DTT ha i paper CS. Servono dati e riferimenti citabili per DTT
e CEFTR. CEFTR: per ora si usa il design point dello script (R0 = 8 m, A = 2.7, B0 = 6 T, 16 TF, ripple 1%); il riferimento pubblicato arriva dall'utente.

## Domande aperte

2. Il confronto out-of-plane del modello globale con ANSYS è stato fatto?
3. Autori, affiliazioni, funding (EUROfusion → disclaimer obbligatorio).
4. Caso mirror di riferimento.

## Compilazione

Serve `iopart.cls` dal template autori IOP (Nuclear Fusion). In alternativa
usare la riga di fallback `article` in `main.tex`.

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```
