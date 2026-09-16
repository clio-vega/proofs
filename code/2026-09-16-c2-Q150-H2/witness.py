"""Explicit two-bead witness.  Geometry (re-derived from eq:2B):
bead A at site 0 moves 0->e ; bead B at site delta moves delta->delta+f.
Order R_f,R_e : A's word has 0 at position delta (B has left), B's word has 0 at e-delta.
Order R_e,R_f : both read 1.  So the element is sigma(pi 0 pi')sigmabar(pi' 0 th')
                                             - sigma(pi 1 pi')sigmabar(pi' 1 th').
For tau = 1_w of rank d, e=d, f=kd: a nonzero element exists exactly when w has period
delta and w_delta = w_{d-delta}  -- the failure condition of Theorem A."""
import sys; sys.path.insert(0,'.')
from independent_maya import commutator, apply_R
from newcounter import tensor_power, Wof
from singleton import crit_star, borders
from itertools import product

def fail_deltas(w):
    n=len(w)
    return [n-p for p in range(n) if w[:p]==w[n-p:] and w[n-p-1]==w[p] and n-p>=1]

def witness_state(w, delta, k):
    d=len(w)+1; e=d; f=k*d
    # sigmabar-word y of length f-1 : blocks w separated by a free bit (set 0)
    y=[]
    for i in range(k):
        y += list(w)
        if i < k-1: y.append(0)
    assert len(y)==f-1
    M=set()
    M.add(0)                                    # bead A
    for i in range(1,delta): 
        if w[i-1]: M.add(i)                     # pi = w_[1,delta-1]
    M.add(delta)                                # bead B
    for i in range(delta+1,d):
        if w[i-1]: M.add(i)                     # pi' = w_[delta+1,d-1]
    # site d must be EMPTY (A's target).  theta' on sites d+1 .. delta+f-1
    for j in range(d-delta+1, f):               # y_j  <-> site delta+j
        if y[j-1]: M.add(delta+j)
    assert delta+f not in M
    return frozenset(M), tuple(y)

def build(d,w,k):
    sig={u:(1 if u==w else 0) for u in product((0,1),repeat=d-1)}
    return Wof(sig,d), Wof(tensor_power(w,d,k),k*d)

print('=== explicit two-bead witnesses for singletons FAILING the border criterion ===')
for d in (4,5,6,7,8,9):
    bad=[x for x in product((0,1),repeat=d-1) if not crit_star(x)]
    shown=0
    for w in bad:
        ds=fail_deltas(w)
        if not ds: print('  !! no failing delta for', w); continue
        delta=ds[0]
        if delta>=d: continue        # delta must be in [1,d-1]
        k=2; W,Wb=build(d,w,k)
        M,y=witness_state(w,delta,k)
        c=commutator({M:1},d,W,k*d,Wb)
        tb=any(len(set(T)^set(M))>=4 for T in c)
        print('  d=%d w=%s delta=%d : commutator %s  (two-bead: %s)'
              % (d,''.join(map(str,w)),delta,'NONZERO '+str(sorted(c.values())) if c else '*** ZERO ***', tb))
        shown+=1
        if shown>=3: break
    sys.stdout.flush()

print()
print('=== same construction on singletons PASSING the criterion: must give zero ===')
for d in (4,5,6,7):
    good=[x for x in product((0,1),repeat=d-1) if crit_star(x)][:3]
    for w in good:
        W,Wb=build(d,w,2)
        n=0
        for delta in range(1,d):
            try: M,y=witness_state(w,delta,2)
            except AssertionError: continue
            if commutator({M:1},d,W,2*d,Wb): n+=1
        print('  d=%d w=%s : %d of %d delta-witnesses nonzero' % (d,''.join(map(str,w)),n,d-1))
    sys.stdout.flush()
