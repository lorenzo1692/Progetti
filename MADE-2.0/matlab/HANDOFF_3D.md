# Handoff aggiornato — modulo 3D, 25 settembre 2026

Leggere prima `docs/TF3D_IMPLEMENTAZIONE.md` e `validation/results/`.

È stato implementato il passo 3D facoltativo del main: geometria dai turn reali, D analitica, fit a tre archi, filamenti planari, Biot–Savart esatto per segmenti a blocchi, campo e carichi campionati, Ampère/ripple, export di curve e centri in mm, grafici. L'iterazione di forma 3D e l'energia approssimata sono opzioni sperimentali disabilitate per default.

**Non confondere il codice consegnato con codice eseguito in MATLAB.** Qui non sono disponibili MATLAB/Octave. Sono stati eseguiti i riferimenti Python indipendenti documentati in `validation/results`; sono consegnati test MATLAB da eseguire. Non è stata rieseguita la regressione FEM 2D. Non sono stati fatti commit o push: il materiale ricevuto è uno zip, senza checkout del repository.

Default cambiati rispetto alla proposta dell'handoff: `shape_mode=0`, `use_arcs=0`, `n_seg=180`. Il fit a tre archi del riferimento Design 7 ha circa 22 mm di scarto e fallisce la tolleranza iniziale di 10 mm, conservata. L'export a tre archi è bloccato se quel criterio fallisce; il percorso analitico resta utilizzabile per lo studio.

Distinzioni essenziali: picco 3D = stima scalare da centro 3D e incremento 2D; quattro sorgenti di quadratura del cavo ≠ cavo raccordato volumetrico; Fz della metà superiore / 2 = media dei carichi dei due tagli, non automaticamente forza della gamba interna; inviluppo radiale del case ≠ solido 3D con wedging e clearances.

Avvio:

```matlab
addpath(genpath(pwd));
validate_tf3d(false);
validate_tf3d(true);
out = run_tf3d_design7();
```

Non usare il risultato per aggiornare il dimensionamento finché non sono chiusi i gate nativi, la convergenza delle forze, il confronto 3D e le ipotesi del case. La modalità iterata richiede anche il test del limite di molte bobine. `out.validated=false` viene mantenuto esplicitamente.

Sono rimasti invariati i sorgenti di scan, campo 2D, FEM 2D e legacy. Restano aperti: riduzione meccanica nello scan; scelta Sm/classificazione; parentesi di B_PHI_TF; progettazione case esterno; confronto di energia/induttanza e carichi con FEM. Non attribuire alla nuova consegna le validazioni dichiarate dal vecchio handoff.


## Aggiornamento: conservazione e salvataggio plot

Tutti i pannelli del FEM surrogato sono conservati: mappa Tresca della sezione e SCL del case, zoom WP con picco jacket, confronto delle tensioni per layer e limiti. La funzione restituisce anche `fig`, senza cambiare le chiamate che non richiedono output. Il main memorizza `mech_fig` e alla fine chiede se salvare tutti i plot aperti, inclusi quelli del FEM e del 3D eseguiti. Rispondere `y`, poi Invio per la cartella proposta oppure digitare un percorso. Si salvano PNG a 300 dpi e FIG modificabili; le figure restano aperte. I file preesistenti non vengono sovrascritti.

Comando manuale per tutti i plot: `ask_save_plots();`
Solo FEM dopo il main: `ask_save_plots(mech_fig);`
Se `mech` esiste ma il plot è stato chiuso: `mech_fig = plot_wp_mech_surrogate(mech,p); ask_save_plots(mech_fig);`

La modifica riguarda grafici e I/O, non i calcoli. Controllo statico eseguito; salvataggio nativo MATLAB non eseguito in questo ambiente.

## Aggiornamento successivo: figura GPS autonoma e uscite EM 3D

Leggere `docs/AGGIORNAMENTO_GPS_EM3D.md`. Il main conserva ora `mech_gps_fig` oltre a `mech_fig`; `plot_wp_gps_section(mech)` permette di rigenerare la sezione dalla soluzione già calcolata. Il plot 3D contiene tutti i percorsi turn, campo su bobina 1 completa per simmetria, profilo radiale e forze. `export_tf3d_results` salva MAT e CSV prima del rendering e indipendentemente dal gate CAD. I limiti fisici di picco, sezione sorgente, forma iterata e induttanza rimangono espliciti. Nessun solver 2D è stato modificato; grafici e test MATLAB non eseguiti qui.
