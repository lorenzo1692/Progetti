import numpy as np, pickle, sys
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
from fem_lin import comp, tresca
tag, femdir, seccsv = sys.argv[1:4]
xs,el=pickle.load(open(f'{tag}_mesh.pkl','rb'))
J=comp(f'{femdir}/TFBM_jacket.txt')
nid=np.array(list(J)); S=np.array([J[n] for n in nid]); XY=np.array([xs[n] for n in nid])
# element-based triangulation of the jacket (quads split), so interpolation never bridges the cable hole
tris=[]
idx={n:i for i,n in enumerate(nid)}
for e,(m,t,r,es,nn) in el.items():
    if m!=2 or len(nn)!=8: continue
    c=[idx.get(n) for n in nn[:4]]
    if None in c: continue
    if c[2]==c[3]: tris.append([c[0],c[1],c[2]])
    else: tris.append([c[0],c[1],c[2]]); tris.append([c[0],c[2],c[3]])
from matplotlib.tri import Triangulation, LinearTriInterpolator
T=Triangulation(XY[:,0],XY[:,1],np.array(tris))
itp=[LinearTriInterpolator(T,S[:,c]) for c in range(4)]
D=np.loadtxt(seccsv,delimiter=',')
res=[]
for r in D:
    a=r[4:6]; b=r[6:8]
    s=np.linspace(0.01,0.99,25); P=a+np.outer(s,b-a)
    V=np.stack([np.asarray(f(P[:,0],P[:,1]).filled(np.nan)) for f in itp],1)
    ok=~np.isnan(V).any(1)
    if ok.sum()<20: res.append((np.nan,np.nan)); continue
    s=s[ok]; V=V[ok]; t=np.linalg.norm(b-a); z=(s-0.5)*t
    m=np.trapezoid(V,s*t,axis=0)/(t*(s[-1]-s[0]))
    bb=6/t**2*np.trapezoid(V*z[:,None],s*t,axis=0)*(1/(s[-1]-s[0]))**0
    res.append((tresca(m),max(tresca(m+bb),tresca(m-bb))))
res=np.array(res)/1e6
sur=D[:,8:10]/1e6; kind=D[:,3]; layer=D[:,10].astype(int)
for kname,kv in (('straight interior',1),('straight end (fillet start)',0),('fillet 45deg',2)):
    m=(kind==kv)&~np.isnan(res[:,0])
    rp=sur[m,0]/res[m,0]; rb=sur[m,1]/res[m,1]
    print(f'{kname:28s} n={m.sum():4d}  Pm sur/FEM {rp.mean():.3f}+-{rp.std():.3f}   PmPb sur/FEM {rb.mean():.3f}+-{rb.std():.3f}   max FEM Pm/PmPb {res[m,0].max():.0f}/{res[m,1].max():.0f}  sur {sur[m,0].max():.0f}/{sur[m,1].max():.0f}')
print('per layer max over all sections  Pm FEM/sur   PmPb FEM/sur')
for k in np.unique(layer):
    m=(layer==k)&~np.isnan(res[:,0])
    print(f'  L{k:<3} {res[m,0].max():5.0f}/{sur[m,0].max():<5.0f} {res[m,1].max():5.0f}/{sur[m,1].max():<5.0f}')
