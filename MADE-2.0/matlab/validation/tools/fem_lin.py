import re, numpy as np, pickle, csv, sys
def comp(fn):
    out={}; on=False
    for ln in open(fn, errors='ignore'):
        if ln.startswith('    NODE     SX'): on=True; continue
        if on:
            if 'MINIMUM VALUES' in ln: break
            p=ln.split()
            if len(p)>=5 and p[0].isdigit():
                v=re.findall(r'-?\d\.\d+E[+-]\d+|0\.0000', ln[8:])
                if len(v)>=4: out[int(p[0])]=[float(x) for x in v[:4]]
    return out
def disp(fn):
    out={}; on=False
    for ln in open(fn, errors='ignore'):
        if 'NODE       UX' in ln: on=True; continue
        if on:
            if 'MAXIMUM ABSOLUTE' in ln: break
            p=ln.split()
            if len(p)>=3 and p[0].isdigit():
                v=re.findall(r'-?\d\.\d+E[+-]\d+|0\.0000', ln[8:])
                if len(v)>=2: out[int(p[0])]=[float(v[0]),float(v[1])]
    return out
def tresca(s):
    sx,sy,sz,txy=s; c=(sx+sy)/2; r=np.hypot((sx-sy)/2,txy); P=[c+r,c-r,sz]; return max(P)-min(P)
def linearize(nodes_xy, S, P0, P1, tol=4e-4):
    """nodes along the segment P0->P1 (thickness line), within tol of the line"""
    d=P1-P0; L=np.linalg.norm(d); u=d/L; nrm=np.array([-u[1],u[0]])
    rel=nodes_xy-P0; s=rel@u; off=rel@nrm
    cand=(np.abs(off)<3e-3)&(s>-1e-6)&(s<L+1e-6)
    if cand.sum()<3: return None
    # pick the node row (constant offset) closest to the line that spans the thickness
    offs=np.unique(np.round(off[cand],5)); offs=offs[np.argsort(np.abs(offs))]
    m=None
    for o in offs:
        mm=cand&(np.abs(off-o)<2e-5)
        if mm.sum()>=3 and s[mm].min()<0.2*L and s[mm].max()>0.8*L: m=mm; break
    if m is None: return None
    ss=s[m]; SS=S[m]; o=np.argsort(ss); ss=ss[o]; SS=SS[o]
    # merge duplicates in s
    us,inv=np.unique(np.round(ss,6),return_inverse=True)
    SS=np.array([SS[inv==i].mean(0) for i in range(len(us))]); ss=us
    if ss[0]>0.2*L or ss[-1]<0.8*L: return None
    mem=np.trapezoid(SS,ss,axis=0)/(ss[-1]-ss[0])
    zc=ss-(ss[0]+ss[-1])/2; t=ss[-1]-ss[0]
    ben=6/t**2*np.trapezoid(SS*zc[:,None],ss,axis=0)
    return mem, ben
