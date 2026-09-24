import re, numpy as np, pickle, collections
def parse(dirn, tag):
    xs={}; el={}
    on=None
    for ln in open(dirn+'/TFBM_global.txt',errors='ignore'):
        if ln.startswith('    NODE        X'): on='n'; continue
        if 'ELEM    MAT    TYP' in ln: on='e'; continue
        if 'LIST NODAL FORCES' in ln: break
        p=ln.split()
        if on=='n' and len(p)==7 and p[0].isdigit():
            xs[int(p[0])]=(float(p[1]),float(p[2]))
        elif on=='e' and len(p)>=8 and all(t.isdigit() for t in p):
            e,mat,typ,rel,esy,sec=map(int,p[:6])
            el[e]=(mat,typ,rel,esy,[int(v) for v in p[6:]])
    pickle.dump((xs,el),open(f'{tag}_mesh.pkl','wb'))
    c=collections.Counter((v[0],v[1],len(v[4])) for v in el.values()); print(tag,len(xs),sorted(c.items()))
if __name__=='__main__':
    import sys
    # usage: python3 parse_fem.py <folder with TFBM_*.txt> <tag>   ->  <tag>_mesh.pkl
    parse(sys.argv[1], sys.argv[2])
