"""Independent Python checks of TF3D mathematics, NOT execution of MATLAB.

Run: python validation/tools/check_tf3d_reference.py
NumPy/SciPy/openpyxl/matplotlib required. Writes results under validation/results.
The geometry and 3D kernel are repeated here deliberately as a reference.
"""
from pathlib import Path
import json, time
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad_vec
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
import openpyxl

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'validation'/'results'

def field(P,A,C,I,block=384):
    P=np.atleast_2d(P);A=np.atleast_2d(A);C=np.atleast_2d(C)
    I=np.broadcast_to(np.asarray(I).reshape(-1),(len(A),))
    D=C-A;length=np.linalg.norm(D,axis=1);E=D/length[:,None]
    B=np.zeros_like(P,dtype=float)
    for i in range(0,len(P),96):
        Q=P[i:i+96];v=np.zeros_like(Q)
        for j in range(0,len(A),block):
            R=Q[:,None,:]-A[None,j:j+block,:];e=E[j:j+block]
            s=np.einsum('msk,sk->ms',R,e);q=s-length[j:j+block]
            cross=np.cross(e[None,:,:],R);rho2=np.sum(cross**2,axis=2)
            near=rho2<=(64*np.finfo(float).eps*max(1,np.max(length[j:j+block])))**2
            if np.any(near&(s>=0)&(q<=0)): raise ValueError('singular filament')
            ra=np.sqrt(rho2+s*s);rb=np.sqrt(rho2+q*q)
            with np.errstate(divide='ignore',invalid='ignore'):
                f=(s/ra-q/rb)/rho2
                stable=length[j:j+block]*(s+q)/(ra*rb*(s*rb+q*ra))
            mask=(s*q>0)&~near;f[mask]=stable[mask];f[near]=0
            v+=1e-7*np.einsum('ms,msk,s->mk',f,cross,I[j:j+block])
        B[i:i+96]=v
    return B

def shape(r1,r2,n=801):
    t=np.linspace(-np.pi/2,np.pi/2,n);k=.5*np.log(r2/r1)
    r=np.sqrt(r1*r2)*np.exp(-k*np.sin(t))
    z=cumulative_trapezoid(-k*np.sin(t)*r,t,initial=0)
    return np.column_stack((r,z))

def loop(U,n=180):
    A=U[-1];s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(U,axis=0),axis=1))]
    ns=max(8,round(n/2*s[-1]/(s[-1]+A[1])));nl=max(2,round(n/2)-ns)
    Q=PchipInterpolator(s,U)(np.linspace(0,s[-1],ns+1));Q[0]=U[0];Q[-1]=A
    line=np.column_stack((np.full(nl+1,A[0]),np.linspace(0,A[1],nl+1)))
    up=np.vstack((line,Q[-2::-1]));low=up[::-1].copy();low[:,1]*=-1
    P=np.vstack((up,low[1:]));D=np.diff(P,axis=0);L=np.linalg.norm(D,axis=1)
    tangent=D/L[:,None];tv=tangent+np.roll(tangent,1,axis=0);tv/=np.linalg.norm(tv,axis=1)[:,None]
    N=np.column_stack((-tv[:,1],tv[:,0]));N=np.vstack((N,N[0]))
    return P,N,len(up)-1

def fit_arcs(U):
    A=U[-1];B=U[0];span=B[0]-A[0]
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(U,axis=0),axis=1))]
    ref=PchipInterpolator(s,U)(np.linspace(0,s[-1],160))
    def build(x):
        if np.any(np.abs(x)>30):return None
        a=np.pi/(1+np.exp(-x[0]));b=a+(np.pi-a)/(1+np.exp(-x[1]))
        al=np.array([0,a,b,np.pi]);mid=span*np.exp(x[2])
        cr=np.cos(al[:-1])-np.cos(al[1:]);sz=np.sin(al[1:])-np.sin(al[:-1])
        M=np.array([cr[[0,2]],sz[[0,2]]])
        if np.linalg.cond(M)>1e8:return None
        ends=np.linalg.solve(M,np.array([span,-A[1]])-mid*np.array([cr[1],sz[1]]))
        R=np.array([ends[0],mid,ends[1]])
        if np.any(R<=span*1e-4) or np.any(R>100*span):return None
        C=[];J=[A]
        for i in range(3):
            C.append(J[-1]+R[i]*np.array([np.cos(al[i]),-np.sin(al[i])]))
            J.append(J[-1]+R[i]*np.array([cr[i],sz[i]]))
        return R,al,np.array(C),np.array(J)
    def dist(P,g):
        R,al,C,J=g;res=np.full(len(P),np.inf)
        for i in range(3):
            v=P-C[i];a=np.mod(np.arctan2(v[:,1],-v[:,0]),2*np.pi)
            mask=(a>=al[i])&(a<=al[i+1]);d=np.minimum(np.linalg.norm(P-J[i],axis=1),np.linalg.norm(P-J[i+1],axis=1))
            d[mask]=np.abs(np.linalg.norm(v[mask],axis=1)-R[i]);res=np.minimum(res,d)
        return res
    def obj(x):
        g=build(x)
        if g is None:return 1e4+np.sum(np.minimum(np.abs(x),100)**2)
        return np.mean((dist(ref,g)/span)**2)
    best=None
    for a in [.5,.9,1.2]:
        for b in [1.9,2.3,2.7]:
            for f in [.4,.9,1.5]:
                res=minimize(obj,[np.log(a/(np.pi-a)),np.log((b-a)/(np.pi-b)),np.log(f)],method='Nelder-Mead',options={'maxiter':500,'maxfev':1200,'xatol':1e-8,'fatol':1e-12})
                if best is None or res.fun<best.fun:best=res
    g=build(best.x);R,al,C,J=g;points=[]
    for i in range(3):
        a=np.linspace(al[i],al[i+1],101)
        q=C[i]+R[i]*np.column_stack((-np.cos(a),np.sin(a)))
        points.extend(q if i==0 else q[1:])
    points=np.array(points);d2=np.full(len(points),np.inf)
    for i in range(len(U)-1):
        v=U[i+1]-U[i];t=np.clip((points-U[i])@v/(v@v),0,1)
        d2=np.minimum(d2,np.linalg.norm(points-U[i]-t[:,None]*v,axis=1))
    err=max(dist(ref,g).max(),d2.max())
    assert np.max(np.abs(J[-1]-B))<1e-10
    return points[::-1],R,C,al,err

def design():
    w=openpyxl.load_workbook(ROOT/'input'/'WP_TF_input_template.xlsx',data_only=True)
    p={r[1]:r[2] for r in w['Input'].iter_rows(min_row=2,values_only=True) if r[1]}
    h=np.repeat([.034167368089179834,.026104457386404455,.024226622311504957],[4,3,6])
    cw=.04300426499905924;jt=.003;ri=1.259425287356322
    Re=ri-p['dr_plasma_side']-p['GoundIns'];x=[];y=[];layers=[];width=[];height=[]
    for k,ch in enumerate(h):
        x.extend((np.arange(8)-3.5)*cw);y.extend([Re-ch/2]*8);layers.extend([k]*8)
        width.extend([cw-2*jt-2*p['turn_insulation_nominal']*p['Increm']]*8)
        height.extend([ch-2*jt-2*p['turn_insulation_nominal']*p['Increm']]*8)
        Re-=ch+p['INS_grades']
    x=np.array(x);y=np.array(y);rc=np.mean(y);r2=(p['R0']+p['R0']/p['A'])*p['ripple']**(-1/p['n_TF'])+ri-rc
    return p,dict(x=x,y=y,u=rc-y,w=np.array(width),h=np.array(height),layer=np.array(layers),r1=rc,r2=r2,I=64628,NI=104*64628)

def sources(P,N,d,nTF,ng=2):
    A=[];C=[];I=[]
    for c in range(nTF):
        if c==0 or ng==0:
            u=np.concatenate([d['u']+a*d['h']/(2*np.sqrt(3)) for a in [-1,1] for b in [-1,1]])
            v=np.concatenate([d['x']+b*d['w']/(2*np.sqrt(3)) for a in [-1,1] for b in [-1,1]])
            currents=np.full(len(u),d['I']/4)
        else:
            u=[];v=[];currents=[]
            for k in np.unique(d['layer']):
                ix=d['layer']==k;uk=np.mean(d['u'][ix])
                if ng==1:
                    u.append(uk);v.append(0);currents.append(np.sum(ix)*d['I'])
                else:
                    spread=np.sqrt(np.mean(d['x'][ix]**2)+np.mean(d['w'][ix]**2)/12)
                    u.extend([uk]*2);v.extend([-spread,spread]);currents.extend([np.sum(ix)*d['I']/2]*2)
        ph=np.pi/2+c*2*np.pi/nTF;er=np.array([np.cos(ph),np.sin(ph),0]);et=np.array([-np.sin(ph),np.cos(ph),0])
        for uj,vj,ij in zip(u,v,currents):
            q=P+uj*N;X=q[:,0,None]*er+vj*et+np.column_stack((np.zeros((len(q),2)),q[:,1]))
            A.append(X[:-1]);C.append(X[1:]);I.append(np.full(len(X)-1,ij))
    return np.vstack(A),np.vstack(C),np.concatenate(I)

def main():
    start=time.time();OUT.mkdir(exist_ok=True,parents=True);result={'scope':'Independent Python mathematical reference; MATLAB code NOT executed'}
    A=np.array([[0.,0.,-1.]]);C=np.array([[0.,0.,1.]])
    b=field([[1.,0.,0.]],A,C,1)[0];expected=2e-7/np.sqrt(2)
    assert np.allclose(b,[0,expected,0],rtol=1e-13,atol=1e-20)
    assert np.allclose(field([[1.,0.,0.]],C,A,1),-b,rtol=1e-13,atol=1e-20)
    assert np.all(field([[0.,0.,2.]],A,C,1)==0)
    try:field([[0.,0.,0.]],A,C,1);raise AssertionError('singularity not rejected')
    except ValueError:pass
    # Check off-axis finite-segment kernel independently by adaptive integration.
    rng=np.random.default_rng(42);errs=[]
    for _ in range(12):
        a=rng.normal(size=3);c=rng.normal(size=3);p=rng.normal(size=3)*2
        dl=c-a
        integ=quad_vec(lambda t:1e-7*np.cross(dl,p-a-t*dl)/np.linalg.norm(p-a-t*dl)**3,0,1,epsabs=1e-16)[0]
        errs.append(np.linalg.norm(field([p],[a],[c],1)[0]-integ)/np.linalg.norm(integ))
    result['finite_segment_max_rel_error_vs_adaptive_quadrature']=max(errs);assert max(errs)<1e-10
    circ=[]
    for n in [64,128,256]:
        t=np.linspace(0,2*np.pi,n+1);Q=np.column_stack((np.cos(t),np.sin(t),np.zeros(n+1)))
        zs=np.array([0,.5,2]);exact=2*np.pi*1e-7/(1+zs**2)**1.5
        bb=field(np.column_stack((np.zeros((3,2)),zs)),Q[:-1],Q[1:],1)[:,2]
        circ.append(float(np.max(np.abs(bb/exact-1))))
    assert circ[2]<circ[1]/3 and circ[1]<circ[0]/3
    result['circular_loop_max_relative_errors_64_128_256']=circ
    p,d=design();U=shape(d['r1'],d['r2']);fit,R,C,angles,err=fit_arcs(U)
    result['design7']={'r1':d['r1'],'r2':d['r2'],'NI':d['NI'],'arc_max_error_m':err,'arc_radii_m':R.tolist()}
    result['arc_fit_10mm_gate_pass']=bool(err<.01)
    # Dense analytic shape with B=C/r has constant curvature*radius.
    t=np.linspace(-np.pi/2,np.pi/2,801);k=.5*np.log(d['r2']/d['r1'])
    r=U[:,0];rp=-k*np.cos(t)*r;zp=-k*np.sin(t)*r
    rpp=k*np.sin(t)*r+k*k*np.cos(t)**2*r
    zpp=-k*np.cos(t)*r+k*k*np.sin(t)*np.cos(t)*r
    curvature=np.abs(rp*zpp-zp*rpp)/(rp*rp+zp*zp)**1.5
    constant=curvature*r
    assert np.max(np.abs(constant*k-1))<1e-12
    result['analytic_curvature_identity_max_error']=float(np.max(np.abs(constant*k-1)))
    mid_values=[];amp=[]
    for n in [90,180,360]:
        P,N,nh=loop(U,n);A,C,I=sources(P,N,d,int(p['n_TF']))
        mid=field(np.column_stack((-d['x'],d['y'],np.zeros(104))),A,C,I)
        mid_values.append(np.linalg.norm(mid,axis=1))
        phi=np.arange(16*int(p['n_TF']))*2*np.pi/(16*p['n_TF']);R0=p['R0']
        b=field(np.column_stack((R0*np.cos(phi),R0*np.sin(phi),np.zeros(len(phi)))),A,C,I)
        bp=-b[:,0]*np.sin(phi)+b[:,1]*np.cos(phi);ref=2e-7*p['n_TF']*d['NI']/R0
        amp.append(float(abs(np.mean(bp)/ref-1)))
        print('reference nseg',n,'sources',len(A),'Ampere err',amp[-1],flush=True)
    assert max(amp)<1e-3
    result['ampere_rel_errors_90_180_360']=amp
    result['centre_field_resolution_max_relative_change_90_to_180']=float(np.max(np.abs(mid_values[0]/mid_values[1]-1)))
    result['centre_field_resolution_max_relative_change_180_to_360']=float(np.max(np.abs(mid_values[1]/mid_values[2]-1)))
    # Compare grouped-neighbour model with all-turn quadrature at z=0.
    P,N,nh=loop(U,180);AA,CC,II=sources(P,N,d,int(p['n_TF']),ng=0)
    mid_full=field(np.column_stack((-d['x'],d['y'],np.zeros(104))),AA,CC,II)
    result['grouped_neighbour_centre_max_relative_error']=float(np.max(np.abs(mid_values[1]/np.linalg.norm(mid_full,axis=1)-1)))
    np.savetxt(OUT/'tf3d_design7_reference_centres.csv',np.column_stack((np.arange(1,105),d['x'],d['y'],*mid_values,np.linalg.norm(mid_full,axis=1))),delimiter=',',header='turn,x2D_m,r_m,Bcentre_90_T,Bcentre_180_T,Bcentre_360_T,Bcentre_full_neighbours_180_T',comments='')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(11,5))
    ax[0].plot(U[:,0],U[:,1],label='Princeton D');ax[0].plot(fit[:,0],fit[:,1],'--',label='Three arcs');ax[0].set(xlabel='R [m]',ylabel='Z [m]',title='Design 7: upper current line');ax[0].axis('equal');ax[0].legend();ax[0].grid(alpha=.25)
    for n,b in zip([90,180,360],mid_values):ax[1].plot(np.arange(1,105),b,label=str(n)+' segments')
    ax[1].set(xlabel='Turn',ylabel='Centre |B| [T]',title='3D centre field at inner mid-plane');ax[1].legend();ax[1].grid(alpha=.25)
    fig.suptitle('Independent Python reference — not MATLAB execution');fig.tight_layout();fig.savefig(OUT/'tf3d_reference.png',dpi=160);plt.close(fig)
    result['elapsed_seconds']=time.time()-start;result['status']='Kernel checks passed; 10 mm arc-fit gate FAILED; MATLAB and FEM validation pending'
    (OUT/'tf3d_python_reference.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
