# Modello a telaio (opzione A): studio per sviluppi futuri

Stato: impostazione (03/10/2026). Non c'è ancora codice. Scopo: un modello
analitico che rappresenti l'effetto telaio del case, cioè quello che il
modello a strati assialsimmetrico non può vedere (`docs/MODELLO_A_STRATI.md`
§10). Oggi la scansione usa al suo posto il surrogato del case (opzione B,
§11 dello stesso documento).

## 1. Cosa deve riprodurre (dal confronto affiancato con il FE)

1. **Pareti laterali tese in direzione radiale** (+30…+60 MPa). Il WP
   spinge sul nose, il nose cede verso l'interno e le pareti lo appendono
   alla piastra lato plasma. Il modello a strati le dà compresse
   (−70…−190 MPa), perché le pareti hanno la stessa deformazione radiale
   del WP.
2. **Ripartizione della forza circonferenziale**: nel FE il nose ne porta di
   più. Quando il nose è sottile la differenza è grande: 25% della forza nel
   nose contro il 18% del modello per s2_4.
3. **Pressione radiale nel WP**. Il WP scarica verso le pareti (attrito,
   scalini). Al centro la pressione è circa il doppio della media sulla
   larghezza, e nei layer profondi cala verso il nose.
4. **Tensione assiale**: più alta nelle pareti (+40%) e nel nose.

## 2. Idea: due campi radiali invece di uno

Nel modello a strati ogni anello ha un solo spostamento radiale u(r),
condiviso da pareti e WP. Nel FE invece le pareti e il WP scorrono uno
rispetto all'altro. Il modello a telaio usa quindi due campi.

- **u_c(r): il case**, cioè nose, pareti laterali e piastra. Le pareti sono
  parallele ai fianchi e i fianchi sono vincolati a restare nel loro piano
  (simmetria ciclica, contatto di wedge senza attrito). Per questo lo
  spostamento radiale del case fissa, a ogni raggio, l'accorciamento della
  corda del settore:

      Δc(r) = 2 u_c(r) tan(π/n_TF)

  Questa è la deformazione che chiude il percorso circonferenziale.
- **u_w(r): il WP**, una colonna caricata dalla forza di Lorentz che
  appoggia sul nose (a r = R_j) e sulle spalle degli scalini, e scorre sulle
  pareti con attrito.

Il percorso circonferenziale a ogni raggio resta quello del modello a
strati, con pareti, isolante e WP in serie lungo la corda. Lo comanda però
u_c, non u_w:

    σ_θ(r) = E_θ,serie(r) · u_c(r)/r

È lo stesso σ_θ per pareti e WP, ed è nel piano di mezzeria che si vede il
WP.

## 3. Equazioni (prima formulazione)

### 3.1 Case

Nel tratto del WP (R_j < r < R_p) il case sono le due pareti, area per
unità di altezza A_w(r) = 2 w(r). La rigidezza radiale viene dalle sole
pareti (non dal WP), quella circonferenziale dalla corda in serie. Equilibrio
radiale di un anello di case, per unità di lunghezza assiale e per settore:

    d(N_r)/dr − (σ_θ chord(r))/r · 1 + τ(r) = 0,     N_r = E_c A_w(r) du_c/dr

dove τ(r) è la forza per unità di raggio che il WP cede alle pareti:

    τ(r) = 2 μ |σ_θ| (attrito) + carichi concentrati sulle spalle dei layer più stretti

Nose (r < R_j) e piastra (r > R_p) restano anelli di Lamé (acciaio
isotropo con l'isolante di wedge), come nel modello a strati, con la
pressione del WP p_w(R_j) applicata in cima al nose.

### 3.2 WP (colonna)

Equilibrio radiale della colonna di larghezza W(r):

    d(p_w W)/dr = f_L(r) W − 2 σ_θ tan(π/n_TF) h-density − τ(r)

- f_L: la forza di Lorentz per unità di volume.
- Il termine in σ_θ è la componente radiale della spinta delle pareti
  inclinate: il termine (σ_r − σ_θ)/r dell'anello.
- Condizioni al contorno: p_w = 0 in cima (lato plasma); la colonna
  appoggia sul nose a R_j.

La deformazione radiale della colonna, ε_r = −p_w/E_r,WP, serve solo per il
jacket. Il carico che il WP trasmette dipende dall'equilibrio, non dalla
sua rigidezza.

### 3.3 Accoppiamenti e incognite

- **Nose.** Il contatto WP–nose a R_j trasferisce p_w(R_j)·W. In più il nose
  si inflette tra le pareti: è una trave curva di luce W + 2 GIT, incastrata
  nelle pareti, che aggiunge al nose al centro la flessione e la membrana
  dell'arco. Questo è il termine che manca al modello a strati quando il
  nose è sottile e il WP largo.
- **Pareti.** Al vault sono collegate al nose e in alto alla piastra; si
  scambiano N_r con il nose e con la piastra, in tensione.
- **Incognite.** u_c(r) per anello come nel modello a strati: due costanti
  per anello più ε_0. Si aggiunge p_w(r), ottenuto in forma chiusa per
  ogni layer dalla ricorsione della colonna (§7 di MODELLO_A_STRATI),
  accoppiato a σ_θ(u_c). Il sistema resta lineare se l'attrito è
  linearizzato (|σ_θ| con il segno noto: compressione).

## 4. Piano di verifica

1. **Stessi controlli del modello a strati**: equilibrio circonferenziale
   e assiale, convergenza con il numero di anelli.
2. **Confronto affiancato** (`validation/compare_layered_vs_fe.m`, da
   estendere) sui 4 design d10, bench, s1_7, s2_4:
   - ripartizione della forza circonferenziale;
   - segno e ordine di grandezza della tensione radiale delle pareti;
   - pressione radiale media del WP per layer: aggiungere al FE la media
     su tutta la larghezza del layer, non solo sulla mezzeria;
   - membrana del nose.
3. **Le 19 run FE di riferimento**, con il criterio usato per il modello
   a strati: case entro ±5% senza correzioni tarate, jacket rms ≤ 10%.
4. **Se passa**, sostituisce nel surrogato del case la correzione
   g(h_WP/t_nose). Se il nose sottile resta fuori, si aggiunge la trave
   curva del nose (§3.3) come termine esplicito.

## 5. Rischi

- La flessione del nose come trave curva tra le pareti è un effetto 2D: in
  un modello radiale entra solo come termine aggiunto, con una luce e un
  incastro da scegliere.
- La pressione del WP non è uniforme sulla larghezza (§1.3): la colonna dà
  la media. Per il jacket serve la distribuzione, ed è lì che il surrogato
  del jacket resta preferibile.
- L'attrito e gli scalini rendono il problema non lineare (stato di
  contatto). Si parte con μ efficace e si verifica la sensibilità.
