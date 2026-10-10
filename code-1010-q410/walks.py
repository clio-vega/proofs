"""HKKO prop:cssyt=path -- K^alpha_{la[m,w]'} = # paths in Z^m from 0 to (la_1..la_m)
staying in { x_1>=x_2>=...>=x_m >= x_1-w }, step i a 0/1 vector with alpha_i ones.
Completely independent of the cylindric Jacobi-Trudi determinant."""
import itertools
from functools import lru_cache

def step_vectors(m, a):
    """0/1 vectors of length m with exactly a ones."""
    out=[]
    for S in itertools.combinations(range(m), a):
        v=[0]*m
        for i in S: v[i]=1
        out.append(tuple(v))
    return out

def count_paths_by_endpoint(m, w, alpha):
    """alpha: tuple of step sizes. Returns dict endpoint -> count."""
    cur={tuple([0]*m):1}
    for a in alpha:
        nxt={}
        for pos,c in cur.items():
            for v in step_vectors(m,a):
                q=tuple(pos[i]+v[i] for i in range(m))
                if all(q[i]>=q[i+1] for i in range(m-1)) and q[m-1]>=q[0]-w:
                    nxt[q]=nxt.get(q,0)+c
        cur=nxt
    return cur

def AB_from_walks(k, w, alpha):
    """A = #tableaux with c^-=+1, B = with c^-=-1, TOT = all of Par(2k,w), content alpha."""
    m=2*k
    ends=count_paths_by_endpoint(m,w,alpha)
    A=B=TOT=0
    for la,c in ends.items():
        TOT+=c
        if all(la[2*i]==la[2*i+1] for i in range(k)): A+=c
        elif la[0]-la[m-1]==w and all(la[2*i+1]==la[2*i+2] for i in range(k-1)): B+=c
    return A,B,TOT
