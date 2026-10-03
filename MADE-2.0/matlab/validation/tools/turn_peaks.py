import numpy as np, pickle, re, sys
from fem_lin import tresca
from fem_lin_generic import cells
tag,femdir,surmat=sys.argv[1:4]
xs,el=pickle.load(open(f'{tag}_mesh.pkl','rb'))
P={}; on=False
for ln in open(f'{femdir}/TFBM_jacket.txt',errors='ignore'):
    if ln.startswith('    NODE     S1'): on=True; continue
    if on:
        if 'MINIMUM VALUES' in ln: break
        v=re.findall(r'-?\d\.\d+E[+-]\d+', ln)
        if len(v)==5: P[int(ln.split()[0])]=float(v[3])
nid=np.array(list(P)); s=np.array([P[n] for n in nid])/1e6; XY=np.array([xs[n] for n in nid])
geo=eval(sys.argv[4]); tins=0.001; JT=float(sys.argv[5])
import scipy.io
m=scipy.io.loadmat(surmat,squeeze_me=True,struct_as_record=False)['out']
sp={(int(t.layer),int(t.col)):(t.peak/1e6,t.peak_xy) for t in m.turn}
rows=[]
for (k,j,x0,x1,y0,y1) in cells(**geo):
    msk=(XY[:,0]>x0)&(XY[:,0]<x1)&(XY[:,1]>y0)&(XY[:,1]<y1)
    i=np.argmax(np.where(msk,s,-1)); rows.append((k,j,s[i],XY[i],sp[(k,j)][0],sp[(k,j)][1]))
print('layer: FEM/sur peak per column')
for k in sorted(set(r[0] for r in rows)):
    rr=[r for r in rows if r[0]==k]
    print(f'L{k:<3}'+' '.join(f'{r[2]:5.0f}/{r[4]:<5.0f}' for r in rr))
for k in (2,4,13):
    r=max([r for r in rows if r[0]==k], key=lambda r:r[2])
    print(f'L{k} FEM max {r[2]:.0f} at {np.round(r[3],4)} col {r[1]} | sur same turn {r[4]:.0f} at {np.round(r[5],4)}')
