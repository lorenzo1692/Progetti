import numpy as np, pickle, sys
from scipy.interpolate import griddata
from fem_lin import comp, tresca
def lin_line(XY, S, P0, P1, n=41):
    s=np.linspace(0,1,n); P=P0+np.outer(s,P1-P0)
    V=np.stack([griddata(XY,S[:,c],P,method='linear') for c in range(4)],1)
    ok=~np.isnan(V).any(1)
    if ok.sum()<n*0.8: return None
    s=s[ok]; V=V[ok]; L=np.linalg.norm(P1-P0); z=(s-0.5)*L
    m=np.trapezoid(V,s*L,axis=0)/(L*(s[-1]-s[0]))
    b=6/L**2*np.trapezoid(V*z[:,None],s*L,axis=0)
    return tresca(m), max(tresca(m+b),tresca(m-b)), m
def scls(Y1, Rin, Ri, h, tw, xcav):
    def xoff(y): return (y*np.sin(h)-tw)/np.cos(h)
    return {'nose x=0':(np.array([0,Y1]),np.array([0,Rin])),
            'nose x=0.10':(np.array([0.10,Y1]),np.array([0.10,np.sqrt(Rin**2-0.01)])),
            'side y=1.05':(np.array([xcav,1.05]),np.array([xoff(1.05),1.05])),
            'side y=0.90':(np.array([xcav,0.90]),np.array([xoff(0.90),0.90])),
            'plate x=0':(np.array([0,Ri-0.02]),np.array([0,Ri])),
            'vault diag':(np.array([xcav,Y1]),np.array([Rin*np.sin(h)-tw*np.cos(h)*0,Rin*np.cos(h)]))}
if __name__=='__main__':
    tag=sys.argv[1]; surcsv=sys.argv[2]; femdir=sys.argv[3]; geo=eval(sys.argv[4])
    xs,el=pickle.load(open(f'{tag}_mesh.pkl','rb'))
    C=comp(f'{femdir}/TFBM_case.txt')
    nid=np.array(list(C)); FS=np.array([C[n] for n in nid]); FX=np.array([xs[n] for n in nid])
    D=np.loadtxt(surcsv,delimiter=','); SX=D[:,1:3]; SS=D[:,3:7]
    fT=np.array([tresca(r) for r in FS])
    print('case nodal SINT max: FEM %.0f at %s | surrogate %.0f at %s'%(fT.max()/1e6,FX[fT.argmax()].round(4),D[:,7].max()/1e6,SX[D[:,7].argmax()].round(4)))
    print('%-14s %9s %9s   %9s %9s'%('SCL','Pm FEM','Pm sur','PmPb FEM','PmPb sur'))
    for name,(a,b) in scls(**geo).items():
        rf=lin_line(FX,FS,a,b); rs=lin_line(SX,SS,a,b)
        if rf is None or rs is None: print(name,'n/a'); continue
        print('%-14s %9.0f %9.0f   %9.0f %9.0f'%(name,rf[0]/1e6,rs[0]/1e6,rf[1]/1e6,rs[1]/1e6))
