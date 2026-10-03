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

Scelte da fare: frame locale (non Frenet, usare il frame finite-build della
bobina), solutore FE proprio in MATLAB vs codice FE esterno pilotato da MADE.

## Codice di riferimento (branch)

| Parte | Branch | Cartella |
|---|---|---|
| TF (versione evoluta), conduttore, modulo 3D, modello globale beam+shell | `claude/code-improvement-3cmng0` | `MADE-2.0/matlab/` (`physics/`, `conductor/`, `coil3d/`, `validation/`) |
| CS e PF | `claude/tf-cs-pfc-folders-ih6f07` | `MADE-2.0/matlab/CS`, `MADE-2.0/matlab/PFC` |
| Vecchi script per macchina (CEFTR 2025, VNS, DEMO…) | `master` | `MADE-2.0/matlab/*.m` |

## Stato dei problemi nel codice (verificato il 2026-10-03 sui branch sopra)

| # | Problema | TF evoluto (`code-improvement`) | CS / PFC | Vecchi script (`master`) |
|---|---|---|---|---|
| 1 | Campo di Ampère `2*pi*RTFi - dr_plasma_side` invece di `2*pi*(RTFi - dr_plasma_side)` (FED 2023 eq. 3 usa `2πR_i`) | **Ancora nel codice** (`physics/compute_operating_params.m:29`). Con `field_model = 1` (default, calibrato) o `2` (discreto) il campo di dimensionamento **non ne dipende**: `k = B_discreto / B_Amp`, quindi `B_Amp` si cancella. Resta nel modello `smeared` (`field_model = 0`), nell'output `B_PHI_0` e nel cavo di riferimento della calibrazione (`1.1*B_Amp`) | non applicabile | **Attivo**: `WP_TF_CEFTR_Design_Design_Point_2025.m`, `WP_TF_VNS_*`, `PLOT_TF_CEFTR/VNS` (sottostima ~3%); gli script DEMO usano `2*pi*RTFi` senza `dr` |
| 2 | Soglia NbTi/Nb3Sn: commento 5 T, codice 6 T | **Risolto**: parametro `cp.B_NbTi_max = 6` in `conductor/cicc_params.m` | soglie diverse dal TF: NbTi < **6.2 T**, HTS > **14.5 T** (TF: 6 / 15 T); anche `d_fili` 0.82 vs 1 mm, VF 0.8 vs 0.7 | `cicc.m` commento sbagliato, `CICC_CEFTR.m` 6 / 6.5 T |
| 3 | Minimo di 4 s su `Tau_discharge` cablato | **Ancora presente** (`search/scan_wp_designs.m:217`, `physics/wp_field_calibration.m:65`) | CS/PF usano `L*Iop/V` | presente |
| 4 | Margine operativo implicito | da dichiarare | da dichiarare | — |

`Tau_discharge1` è l'eq. (9) di FED 2023 (`S_VV` = ammissibile del vacuum vessel, 120 MPa): ora è citata nel paper.

Per il paper il punto 2 conta: se diciamo "modulo conduttore condiviso", oggi
in realtà ci sono tre copie con parametri diversi. Conviene unificarle in
`conductor/` con le soglie come input.

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
e CEFTR (il CEFTR nello script: R0 = 8 m, A = 2.7, B0 = 6 T, 16 TF, ripple 1%).

## Domande aperte

1. CEFTR: qual è il riferimento pubblicato del design magnetico?
2. Il confronto out-of-plane del modello globale con ANSYS è stato fatto?
3. Autori, affiliazioni, funding (EUROfusion → disclaimer obbligatorio).
4. Caso mirror di riferimento.

## Compilazione

Serve `iopart.cls` dal template autori IOP (Nuclear Fusion). In alternativa
usare la riga di fallback `article` in `main.tex`.

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```
