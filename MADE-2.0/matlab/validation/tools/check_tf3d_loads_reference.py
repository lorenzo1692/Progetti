"""Additional Python reference: 2D centres and upper-half Lorentz forces.
Does not execute MATLAB or certify the 3D force/peak model.
"""
import numpy as np,json,time
from check_tf3d_reference import design,shape,loop,sources,field,OUT

def centres2d(d,nTF):
    px=d['x'];py=d['y'];bx=np.zeros(len(px));by=bx.copy()
    def G(u,v):
        with np.errstate(divide='ignore',invalid='ignore'):
            t1=.5*v*np.log(u*u+v*v);t2=u*np.arctan(v/u)
        return np.where(v==0,0,t1)+np.where(u==0,0,t2)
    for x,y,w,h in zip(px,py,d['w'],d['h']):
        r=min(.003,.49*min(w,h));J=d['I']/(w*h-(4-np.pi)*r*r)
        rect=[(x,y,w-2*r,h)]
        for sx in [-1,1]:rect.append((x+sx*(w/2-r/2),y,r,h-2*r))
        yy=np.linspace(0,r,9)
        F=lambda y:.5*(y*np.sqrt(np.maximum(0,r*r-y*y))+r*r*np.arcsin(np.minimum(y/r,1)))
        for k in range(8):
            wk=(F(yy[k+1])-F(yy[k]))/(yy[k+1]-yy[k]);yb=(yy[k+1]+yy[k])/2
            for sx in [-1,1]:
                for sy in [-1,1]:rect.append((x+sx*(w/2-r+wk/2),y+sy*(h/2-r+yb),wk,yy[k+1]-yy[k]))
        for x,y,w,h in rect:
            ua=px-(x-w/2);ub=px-(x+w/2);va=py-(y-h/2);vb=py-(y+h/2)
            by+=2e-7*J*(G(ua,va)-G(ua,vb)-G(ub,va)+G(ub,vb))
            bx-=2e-7*J*(G(va,ua)-G(vb,ua)-G(va,ub)+G(vb,ub))
    for c in range(1,nTF):
        a=c*2*np.pi/nTF;x=d['x']*np.cos(a)-d['y']*np.sin(a);y=d['x']*np.sin(a)+d['y']*np.cos(a)
        dx=px[:,None]-x;dy=py[:,None]-y;dd=dx*dx+dy*dy
        bx+=2e-7*d['I']*np.sum(-dy/dd,axis=1);by+=2e-7*d['I']*np.sum(dx/dd,axis=1)
    return np.sqrt(bx*bx+by*by)

def evaluate(n):
    p,d=design();U=shape(d['r1'],d['r2']);P,N,nh=loop(U,n);A,C,I=sources(P,N,d,int(p['n_TF']))
    points=[];dl=[]
    for u,v in zip(d['u'],d['x']):
        q=P+u*N;Q=np.column_stack((-v*np.ones(len(q)),q))
        points.append((Q[:nh]+Q[1:nh+1])/2);dl.append(np.diff(Q[:nh+1],axis=0))
    points=np.vstack(points);dl=np.vstack(dl);B=field(points,A,C,I)
    F=d['I']*np.sum(np.cross(dl,B),axis=0)
    mid=field(np.column_stack((-d['x'],d['y'],np.zeros(104))),A,C,I)
    twod=centres2d(d,int(p['n_TF']));three=np.linalg.norm(mid,axis=1)
    bmag=np.linalg.norm(B,axis=1)
    result=dict(nseg=len(P)-1,n_half=nh,half_force_N=F.tolist(),mean_cut_force_N=float(F[2]/2),centre_max_T=float(bmag.max()),midplane_2d_max_T=float(twod.max()),midplane_3d_max_T=float(three.max()),midplane_max_relative_difference=float(np.max(np.abs(three/twod-1))))
    np.savetxt(OUT/f'tf3d_upper_reference_{n}.csv',np.column_stack((points,B)),delimiter=',',header='X_m,Y_m,Z_m,BX_T,BY_T,BZ_T',comments='')
    return result
if __name__=='__main__':
    start=time.time();a=evaluate(90);print(a,flush=True);b=evaluate(180);print(b,flush=True)
    result=dict(scope='Python reference only, not MATLAB execution',resolutions=[a,b],mean_force_relative_change=abs(a['mean_cut_force_N']/b['mean_cut_force_N']-1),elapsed_seconds=time.time()-start)
    result['force_5percent_resolution_gate_pass']=result['mean_force_relative_change']<.05
    (OUT/'tf3d_loads_python_reference.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
