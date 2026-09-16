"""Side-effect-free helpers for the Q150-H2 session."""
import sys; sys.path.insert(0,'.')
from independent_maya import commutator, apply_R, alpha, L, inM, word
from itertools import product

def borders(w):
    n=len(w); return [p for p in range(n) if w[:p]==w[n-p:]]
def crit_star(w):
    n=len(w); return all(w[p]!=w[n-1-p] for p in borders(w))
def crit_starstar(w):
    n=len(w); d=n+1
    return all(w[dd-1]!=w[d-dd-1] for dd in range(1,d))
def fail_deltas(w):
    """delta in [1,d-1] witnessing failure of crit_star: period delta and w_delta=w_{d-delta}."""
    n=len(w)
    return [n-p for p in range(n) if w[:p]==w[n-p:] and w[n-p-1]==w[p] and n-p>=1]

def tensor_power(w,d,k):
    m=k*d; S={}
    for y in product((0,1),repeat=m-1):
        v=1
        for i in range(k):
            if y[i*d:i*d+d-1]!=w: v=0; break
        S[y]=v
    return S
def Wof(sig): return {u:((-1)**sum(u))*v for u,v in sig.items()}
def build(d,w,k):
    sig={u:(1 if u==w else 0) for u in product((0,1),repeat=d-1)}
    return Wof(sig), Wof(tensor_power(w,d,k))

def witness_state(w,delta,k):
    """Initial Maya set for the two-bead element at delta:
       bead A at 0 moves 0->d, bead B at delta moves delta->delta+kd."""
    d=len(w)+1; f=k*d
    y=[]
    for i in range(k):
        y+=list(w)
        if i<k-1: y.append(0)
    M={0, delta}
    for i in range(1,delta):
        if w[i-1]: M.add(i)
    for i in range(delta+1,d):
        if w[i-1]: M.add(i)
    for j in range(d-delta+1,f):
        if y[j-1]: M.add(delta+j)
    assert d not in M and delta+f not in M
    return frozenset(M), tuple(y)
