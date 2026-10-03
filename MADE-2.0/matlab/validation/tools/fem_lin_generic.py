import numpy as np, pickle, sys
from fem_lin import *
def cells(Ri_, dps, GIT, INS, n_turns, Cond_w, Cond_h):
    out=[]; Re=Ri_-dps-GIT
    for k,(nt,w,h) in enumerate(zip(n_turns,Cond_w,Cond_h)):
        Ri=Re-h
        for j in range(nt):
            x0=-nt*w/2+j*w; out.append((k+1,j+1,x0,x0+w,Ri,Re))
        Re=Ri-INS
    return out
def run(tag, femdir, geo, JT, tins, outpkl, fr=(0.25,0.5,0.75)):
    xs,el=pickle.load(open(f'{tag}_mesh.pkl','rb'))
    J=comp(f'{femdir}/TFBM_jacket.txt')
    nid=np.array(list(J)); S=np.array([J[n] for n in nid]); XY=np.array([xs[n] for n in nid])
    rows=[]
    for (k,j,cx0,cx1,cy0,cy1) in cells(**geo):
        o=tins+JT[k-1]; x0,x1,y0,y1=cx0+o,cx1-o,cy0+o,cy1-o; t=JT[k-1]
        for fi,f in enumerate(fr):
            yc=y0+f*(y1-y0); xc=x0+f*(x1-x0)
            secs={'left':(np.array([x0-t,yc]),np.array([x0,yc])),'right':(np.array([x1,yc]),np.array([x1+t,yc])),
                  'bottom':(np.array([xc,y0-t]),np.array([xc,y0])),'top':(np.array([xc,y1]),np.array([xc,y1+t]))}
            for w,(a,b) in secs.items():
                r=linearize(XY,S,a,b)
                if r is None: continue
                m,bb=r
                rows.append((k,j,w,fi+2,tresca(m),max(tresca(m+bb),tresca(m-bb)),*m))
    pickle.dump(rows,open(outpkl,'wb')); print('FEM sections',len(rows))
if __name__=='__main__':
    # example (first benchmark):
    run('bench','fem_run',dict(Ri_=1.18,dps=0.02,GIT=0.005,INS=0.0005,n_turns=[10]*6+[8]*5,Cond_w=[0.046]*11,Cond_h=[0.0326]*6+[0.0228]*5),
        [0.0035]*11, 0.001, 'femB_quarter.pkl')
