# Sezione GPS e completamento delle uscite EM 3D

La figura autonoma `plot_wp_gps_section` usa direttamente `mech.mesh` e `mech.sol.U` del FEM 2D generalized plane strain. Mostra mesh/materiali, Tresca, deformata e spostamento radiale del WP. Riporta la deformazione assiale globale e lo stato dei controlli numerici. Il fattore di amplificazione della deformata è dichiarato, mentre la scala colore degli spostamenti conserva i millimetri fisici. I materiali mantengono gli ID del solver: 1 LTS, 2 jacket, 3 isolamento, 4 case, 5 wedge insulation, 6 HTS, 7 filler.

Il main ora riceve `[mech_fig, mech_gps_fig]` da `plot_wp_mech_surrogate`: la nuova figura GPS è generata insieme a quella precedente con le SCL e i confronti per layer. La domanda finale di salvataggio include entrambe. Le chiamate precedenti senza output, o con un solo output, rimangono compatibili.

Se il risultato `mech` è già nel workspace, senza ripetere il calcolo:

```matlab
mech_gps_fig = plot_wp_gps_section(mech);
ask_save_plots(mech_gps_fig);
```

Per rigenerare tutte le figure meccaniche:

```matlab
[mech_fig, mech_gps_fig] = plot_wp_mech_surrogate(mech,p);
ask_save_plots([mech_fig, mech_gps_fig]);
```

## Che cosa mancava nella presentazione EM 3D

La precedente figura mostrava solo una linea baricentrica per bobina e i campioni di campo sulla metà superiore della bobina 1. Questa rappresentazione non rendeva visibile il pacchetto completo, benché le sorgenti magnetiche fossero circuiti chiusi. Mancavano inoltre il profilo radiale del campo, i grafici dei carichi e l'export tabellare dei campi vettoriali.

Ora sono presenti:

- percorsi chiusi di ogni turn di tutte le bobine, con offset toroidali planari e senza crossover/terminali;
- campo centrale lungo tutto il percorso della bobina 1: metà superiore calcolata, metà inferiore ricostruita per simmetria;
- curve lungo l'intera ascissa curvilinea, confronto centrale 2D/3D, ripple e profilo radiale Bphi nel volume del plasma;
- vettori di forza per intervallo di integrazione e forza netta di centraggio per turn;
- massimi campionati per turn, posizione dei massimi e distinzione tra campo centrale calcolato e picco stimato;
- matrice di induttanza equivalente e relativa figura quando l'opzione di energia viene esplicitamente abilitata;
- export separato dei risultati numerici, indipendente dall'accettazione del CAD e precedente ai grafici.

`export_tf3d_results` scrive MAT, CSV di campo/forza sulla bobina 1 completa, riepilogo turn, ripple, profilo radiale e stato dei controlli. I CSV riportano quali righe provengono dalla simmetria. Le forze sono espresse in N per intervallo rappresentato, non N/m. L'annullamento della risultante verticale dell'intera bobina deriva dalla simmetria imposta e **non è un controllo indipendente di equilibrio**. Gli inviluppi del case nei grafici ora usano lo stesso nose esterno variabile dell'export.

## Quanto è completo il modello fisico

È stato completato questo insieme di uscite e rappresentazioni; **non è ancora un modello EM 3D validato in tutti i suoi aspetti**. Restano gli stessi limiti fisici: quadratura 2×2 della sezione del cavo, raggruppamento opzionale delle bobine vicine, picco stimato con incremento scalare 2D, forma iterata sperimentale e induttanza da linea baricentrica/GMD. Non sono stati convertiti in risultati esatti, né abilitati automaticamente. Il case visualizzato non è un solido FEM 3D e non si simulano le connessioni fra turn.

Il riferimento Python verifica la parità del campo confrontando direttamente Biot–Savart sopra e sotto il piano equatoriale: errore relativo circa 1,54e-16 su tre punti di controllo. I dati sono in `validation/results/tf3d_completeness_reference.json`. Il test nativo MATLAB include ora questo confronto e il conteggio delle righe dell'export dell'intera bobina.

**MATLAB/Octave non sono disponibili nell'ambiente di consegna:** il rendering dei nuovi grafici MATLAB e i test nativi sono da eseguire. I sorgenti del solver GPS, del campo 2D e dello scan non sono stati modificati; non è stata rieseguita la regressione meccanica. Per una verifica nativa iniziare con `validate_tf3d(false)` e poi `validate_tf3d(true)`, mantenendo anche il benchmark FEM 2D come gate nel proprio ambiente.
