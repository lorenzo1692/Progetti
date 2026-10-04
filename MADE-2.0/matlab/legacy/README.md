# legacy/ — codice non usato dalla pipeline

Qui ci sono gli script e le funzioni della versione precedente del tool,
conservati senza modifiche come riferimento. **Nessun file della pipeline
attuale li chiama** (verificato con un'analisi delle dipendenze a partire
da `main_WP_TF_design.m`, dagli script di `validation/` e da
`postprocess/plot_wp_section_manual_example.m`). I punti di ingresso
tolgono questa cartella dal path (`rmpath(genpath(... 'legacy'))`), così
nessuna funzione vecchia può oscurare quelle nuove.

| Cartella | Contenuto | Sostituito da |
|---|---|---|
| `scripts/` | Script monolitici di progetto per macchina (`WP_TF_*_Design_*.m`, `WP_TF_Esplorazione.m`) e di plot (`PLOT_TF_*.m`, `Plot_WP_TF.m`), con i parametri scritti nel codice e `shape_cable` impostato a mano (201, oppure le stringhe `'RIS_'`/`'Rect'`) | `main_WP_TF_design.m` + file di input Excel; `postprocess/plot_wp_section*.m` |
| `conductor/` | Versioni per macchina del dimensionamento CICC (`CICC_CEFTR.m`, `CICC_DEMO.m`), leggi Ic separate per strand (`Ic_Nb3Sn_DTT_TF_KAT.m`, `Ic_Nb3Sn_ITER.m`, `Ic_Nb3Sn_WST.m`, `Ic_NbTi_TF.m`), `Jc_REBCO.m`, `Ths_TF_Nb3Sn.m`, `scf.m` | `conductor/cicc.m`, `Ic_Nb3Sn.m` (le tre leggi Nb3Sn in una, scelta con il tipo di strand), `Ic_NbTi.m`, `Ic_sst33.m` + `jc_ybco.m`, `heat_balance_cicc_ode.m`; lo SCF è calcolato in `search/scan_wp_designs.m` |
| `emag_axisym/` | Utilità del vecchio calcolo elettromagnetico assialsimmetrico: campo di spire circolari (`xbr.m`, `xbz.m`), induttanze (`xlsheet.m`), `ColorQuiver.m`; `xbr.m`, `xbz.m`, `xlm.m` sono passati in `shared/emag_axisym/` (usati da CS e PFC) | `physics/wp_field_at_points.m` (2D), `coil3d/biot_savart_segments.m` (3D) |
| `coil3d/` | Flusso ANSYS per la shape: `MAIN.dat` (sorgenti SOURC36 lungo la shape, BIOT, campo sulla linea media), `TF_shape_opt_BF_Ncoils.m` e `bendingfree_opt.m` (shape bending-free e fit a 3 archi), `emag_calculation_TFC*.m` | `coil3d/tf_bending_free_shape.m`, `tf3d_iterate_shape.m`, `tf_three_arc_fit.m`, `biot_savart_segments.m` |

## Eseguire un file legacy

Gli script legacy usano ancora alcune funzioni attive (per esempio `cicc`
e `heat_balance_cicc_ode`). Per eseguirne uno, aggiungere al path sia la
cartella principale sia `legacy/`:

```matlab
addpath(genpath('<...>/MADE-2.0/matlab'));   % include legacy/
run('<...>/MADE-2.0/matlab/legacy/scripts/WP_TF_VNS_Design_Point_2026.m')
```

Non modificare questi file: le correzioni vanno fatte nella pipeline attuale.
