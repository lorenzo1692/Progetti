# Passaggio 3D a valle di MADE — 25 settembre 2026

Per la figura GPS autonoma e le successive estensioni delle uscite EM leggere anche `AGGIORNAMENTO_GPS_EM3D.md`.

## Stato della consegna

È stato aggiunto un modulo MATLAB di geometria e magnetostatica 3D, collegato alla soluzione selezionata dal main. È una versione di ricerca da eseguire e verificare in MATLAB: **in questo ambiente non sono disponibili né MATLAB né Octave**. I risultati numerici allegati provengono da un'implementazione Python indipendente delle formule, non dall'esecuzione dei file `.m`. Non è stata eseguita una validazione FEM 3D.

Le funzioni di dimensionamento, campo 2D, meccanica 2D e scan preesistenti sono rimaste identiche ai file dello zip ricevuto. Il main contiene un nuovo passo facoltativo e il template Excel contiene 15 nuovi parametri. Non è stato modificato il denominatore segnalato in `compute_operating_params`, né sono stati reinterpretati gli allowables meccanici. La validazione FEM 2D riportata nel precedente handoff resta un risultato ereditato, non rieseguito qui.

## Uso

Dalla cartella `matlab`:

```matlab
addpath(genpath(pwd));
report = validate_tf3d(false);    % kernel, geometria, archi e Ampere
report = validate_tf3d(true);     % aggiunge caso 7 a due risoluzioni ed export
out = run_tf3d_design7();         % esempio completo, grafici ed export
```

Per una soluzione selezionata o caricata con `load_design_point`:

```matlab
out = tf3d_from_design(row, p);
plot_tf3d(out, p);
files = export_tf3d_shape(out, fullfile(pwd, 'TF3D_output'), 'selected_design');
```

Il main propone il nuovo passo dopo la verifica meccanica facoltativa; l'export ANSYS 2D diventa il passo 7. La funzione di calcolo non scrive file autonomamente: il flag `tf3d_export` viene usato dal main. L'esempio `run_tf3d_design7` esporta esplicitamente.

Gli override in una chiamata prevalgono sull'Excel e omettono il prefisso `tf3d_`:

```matlab
out = tf3d_from_design(row, p, struct('n_seg',360,'n_eval',180, ...
    'neighbor_gauss',0,'use_arcs',0,'shape_mode',0,'energy',0));
```

Questa chiamata è più costosa: conserva la quadratura di tutti i turn anche nelle bobine vicine. Il modulo è a valle dello scan e non è progettato per centomila candidati.

## Geometria e convenzioni

La sezione resta quella di MADE: x toroidale, y radiale; il layer 1 è lato plasma. I centri dei turn sono ricostruiti dalle dimensioni e dal numero di conduttori di ciascun layer. A correnti uguali il centroide della corrente è la media dei raggi dei turn, non il centro del rettangolo che contiene il WP.

Si usa `r1 = mean(y_turn)`, cioè il raggio della **soluzione selezionata**. Il raggio macchina `g.RTFi` viene conservato come metadato, senza sostituirlo a `row.Ri_`. Il raggio esterno della linea di corrente è `r2 = g.RTFo + (row.Ri_ - r1)`. Il WP è quindi specchiato sulla gamba esterna, come ipotesi geometrica del piano di handoff.

Ogni turn è un circuito chiuso indipendente. Sono omessi salti di layer, crossover, terminali e raccordi fra bobine. Nel riferimento globale la bobina 1 giace nel piano YZ; la sua normale toroidale è −X. Pertanto `X = -x_turn` sulla gamba interna. Per le altre bobine gli offset toroidali seguono la normale al piano, senza deformare la bobina con una rotazione dipendente dal raggio.

Il raggio del case sul piano centrale interno è `Rk_/cos(pi/n_TF)`. La proiezione esterna usa per default un nose pari a metà di quello interno; la piastra lato plasma mantiene lo spessore interno. Si tratta di parametri e ipotesi geometriche, non di dimensionamento verificato del case esterno. L'export fornisce inviluppi radiali e linee, non solidi wedged, né verifica l'interferenza fra case di bobine adiacenti.

## Forma libera e archi

Il default è `shape_mode=0`: D analitica nel modello di campo efficace proporzionale a 1/r. La curva viene ricostruita sulla linea di corrente. La gamba interna rettilinea completa il circuito; la condizione di tensione costante riguarda la parte curva libera, non elimina le reazioni della gamba vincolata.

La modalità `shape_mode=1` implementa un punto fisso sperimentale: campiona il campo 3D medio efficace sui turn, aggiorna la forma mediante l'integrale radiale e applica un rilassamento. Il residuo è quello della proposta non rilassata. Un arresto per numero massimo di iterazioni non viene dichiarato convergenza. **Questa modalità non è stata eseguita né validata nella consegna**, incluso il limite di molte bobine previsto dall'handoff. Resta disabilitata per default.

Il fit a tre archi è deterministico. Tre intervalli di rotazione della tangente, con raggi positivi, garantiscono continuità di posizione e tangente; due raggi sono ottenuti dalle condizioni esatte sui due estremi. La ricerca usa `fminsearch` con una griglia iniziale fissa e un numero massimo di valutazioni. Lo scarto è misurato in entrambi i sensi con campionamento di curve e segmenti.

Sul riferimento Python del Design 7, lo scarto massimo è circa **22,01 mm**. Il criterio iniziale di 10 mm non è stato allargato. Il default è quindi `use_arcs=0`: il calcolo usa la curva libera. Se si seleziona `use_arcs=1`, questo fit non supera il controllo e l'export CAD viene bloccato. La scelta di una tolleranza costruttiva diversa richiede una decisione progettuale, oppure un fit/una rappresentazione più ricca. I tre archi vengono comunque riportati come confronto.

## Campo 3D e auto-campo

`biot_savart_segments` integra esattamente il contributo di ciascun segmento rettilineo. Usa proiezione longitudinale, distanza dalla retta e una forma razionalizzata per i punti esterni al segmento, evitando cancellazioni numeriche. Il calcolo usa blocchi di punti e sorgenti: non alloca la matrice completa punto × segmento.

Un punto appartenente a un segmento sorgente è singolare e genera un errore. Non viene sostituito con campo nullo. I punti sulla retta ma esterni al segmento hanno invece contributo nullo, come previsto dalla geometria del prodotto vettoriale.

Per valutare i centri dei cavi senza la singolarità del filamento centrale, la bobina 1 usa quattro filamenti di quadratura per ogni cavo (2 × 2, sezione rettangolare equivalente). Questo rende il campo centrale definito e conserva un contributo della curvatura del proprio circuito. **Non equivale all'integrazione volumetrica del cavo raccordato**: l'errore della quadratura della sezione deve essere studiato, soprattutto sulle curve e per i carichi.

Per le bobine vicine, il default usa due sorgenti per layer che conservano corrente, centroide e secondo momento toroidale della distribuzione rappresentata. `neighbor_gauss=0` conserva tutti i turn con quadratura. Il confronto Python sul piano equatoriale del Design 7 dà uno scarto massimo relativo del campo centrale di circa **0,527%** tra i due modelli, con lo scarto relativo maggiore sui turn a basso campo. Non è una dimostrazione dell'errore massimo su tutta la bobina.

Il picco riportato è una **stima scalare**:

\[
B_{peak,est}(s)=|B_{centre,3D}(s)|+
[B_{peak,2D}-|B_{centre,2D}|].
\]

Non è una sovrapposizione vettoriale, né una ricerca 3D del massimo sulla sezione. L'incremento 2D viene trasportato lungo il percorso senza una validazione sulle curve; questa stima non aggiorna `B_grade` o il dimensionamento del superconduttore. I massimi lungo il percorso sono massimi dei campioni, non massimi continui garantiti.

## Carichi e forza longitudinale

I carichi vengono calcolati con `I dl × B` sui centri dei turn della metà superiore. Gli incrementi di lunghezza sono quelli dei percorsi offset; con campionamento ridotto ogni punto rappresenta una porzione della linea. `n_eval` deve quindi essere aumentato insieme a `n_seg` per verificare le forze, non solo il campo.

Per il corpo libero della metà superiore:

\[
F_{z,upper}=T_{inner}+T_{outer}.
\]

Di conseguenza `Fz_upper/2` è la **media** delle due forze di taglio assiali. L'uguaglianza con la tensione della gamba interna richiede un'ulteriore ipotesi di tensione costante. Il valore è confrontato con la formula MADE, ma non viene passato automaticamente al solver meccanico 2D. La forma a tre archi e il pacchetto finito possono sviluppare bending.

Il campo medio toroidale su un cerchio che concatena tutti i turn interni viene confrontato con Ampère usando gli ampere-spire effettivi, non quelli richiesti prima dell'arrotondamento. La legge integrale è esatta, mentre la media campionata ha una tolleranza numerica. Il ripple è campionato su un settore al bordo esterno del plasma, sia per Bphi sia per |B|.

## Energia

`energy=1` attiva una stima opzionale di Neumann con una sola linea baricentrica per bobina, corrente NI e un kernel regolarizzato con GMD approssimata del WP rettangolare. La matrice riguarda questi circuiti equivalenti, non i singoli turn. L'energia contiene il fattore N² e si riporta `L_eq = 2W/(n_TF Iop²)`.

Questa stima non è stata validata, resta disattivata per default e non modifica scarica o hot spot. Una discretizzazione convergente del kernel regolarizzato non dimostrerebbe comunque la correttezza dell'approssimazione di sezione.

## Controlli e limiti della verifica eseguita

`validation/results/tf3d_python_reference.json` contiene i risultati riproducibili del controllo indipendente, con CSV e grafico. I controlli comprendono segmento finito contro quadratura adattativa, orientamento, singolarità, spira circolare, curvatura della D analitica, tangenza degli archi, media di Ampère e sensibilità del campo centrale al raggruppamento delle bobine vicine.

Il controllo sulla spira circolare dà errori massimi relativi di 0,11246%, 0,02811% e 0,00703% con 64, 128 e 256 segmenti. Il kernel di segmento concorda con la quadratura adattativa entro circa 7,4 × 10⁻¹⁶ nei dodici punti del test. Questi numeri non sono l'accuratezza complessiva del modello TF.

Un secondo riferimento Python ha integrato i carichi su tutti i segmenti della metà superiore. Con 90 e 180 segmenti per circuito, la media delle due forze assiali di taglio è rispettivamente **44,652 e 44,877 MN**, con variazione dello **0,502%**. Il massimo campo centrale campionato sull'intera metà superiore passa da 13,814 a 13,869 T. Il campo centrale sul piano equatoriale differisce dal modello 2D raccordato di al massimo **1,742%** sui turn: il massimo del campo centrale è 13,733 T nel 3D e 13,737 T nel 2D. Sono confronti del **campo al centro**, non del picco sulla superficie del cavo. Queste prove usano la forma analitica e due punti di quadratura per layer delle bobine vicine; non separano ancora l'errore della quadratura del cavo da quello della geometria 3D. I risultati sono in `tf3d_loads_python_reference.json`.

Il test nativo `validate_tf3d` è consegnato ma non eseguito. Esso separa la correttezza del codice dal criterio geometrico sugli archi, che resta esplicitamente fallito sul Design 7. Anche se i controlli numerici passano, `out.validated` resta `false`: manca la validazione 3D. `out.status` non contiene un verdetto di ammissibilità strutturale.

La 2D FE regression non è stata rieseguita: non è stato modificato alcun suo sorgente, ma questo non sostituisce il gate nativo richiesto dall'handoff prima di integrare la consegna nel repository.

## File principali

| File | Funzione |
|---|---|
| `coil3d/tf3d_from_design.m` | Driver a valle della soluzione MADE |
| `coil3d/tf3d_design_geometry.m` | Centri reali, raggi e inviluppi |
| `coil3d/tf_bending_free_shape.m` | Forma analitica o da campo efficace |
| `coil3d/tf3d_iterate_shape.m` | Iterazione 3D sperimentale |
| `coil3d/tf_three_arc_fit.m` | Fit deterministico e scarto |
| `coil3d/tf3d_coil_filaments.m` | Sorgenti planari con sezione finita approssimata |
| `coil3d/biot_savart_segments.m` | Kernel analitico a blocchi |
| `coil3d/tf3d_evaluate_turns.m` | Campo centrale e carichi campionati |
| `coil3d/tf3d_toroidal_checks.m` | Ampère e ripple |
| `coil3d/export_tf3d_shape.m` | CSV CAD in mm, centri turn, risultati MAT |
| `postprocess/plot_tf3d.m` | Quattro viste di geometria e campo |
| `validation/validate_tf3d.m` | Test nativi da eseguire |
| `validation/tools/check_tf3d_reference.py` | Controllo matematico Python indipendente |
| `validation/tools/check_tf3d_loads_reference.py` | Confronto di centri 2D e forze 3D |

## Passi successivi

Prima eseguire i test nativi e il benchmark meccanico ereditato nel proprio ambiente MATLAB/Octave. Poi confrontare geometria, campo vettoriale e risultanti contro ASPIRE/ANSYS 3D su una geometria fissa; variare discretizzazione del percorso, quadratura della sezione e rappresentazione delle bobine vicine separatamente. Solo dopo usare il campo per un'iterazione di forma o aggiornare i carichi del modello meccanico.

La riduzione meccanica a componenti condensati da usare **dentro lo scan** resta un'attività distinta e aperta: non è stata sostituita da questo modulo 3D.

Riferimento consultato per la convenzione di Biot–Savart e la verifica del segmento: [OpenStax, University Physics 2, §12.1](https://openstax.org/books/university-physics-volume-2/pages/12-1-the-biot-savart-law). Il codice, i criteri e i risultati numerici di questa consegna sono originali; il piano progettuale di partenza è `HANDOFF.md` fornito dall'utente.
