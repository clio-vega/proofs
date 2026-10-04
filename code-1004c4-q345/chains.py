"""Side A: Samuel's T_{w/u} (MO 313951).

V = span{ t_{ab} : a<b } / ( t_{ab} + t_{bc} = t_{ac} ).
Writing e_i := t_{i,i+1} one gets t_{ab} = e_a + ... + e_{b-1}, i.e. V is FREE on
the e_i and t_{ab} |-> the positive root alpha_a + ... + alpha_{b-1}.  So an element
of V^{tensor m} is a dict  (p_1,...,p_m) -> integer,  p_i in {1,...,n-1}.

T_{w/u} = sum over factorisations u * t_1 * ... * t_m = w with the prefix lengths
increasing by 1 at every step, of t_1 (x) ... (x) t_m.
"""
from collections import Counter
from schub import length, rmul_t, refl_length

def covers(u, n, lenfn=length):
    """all (a,b,ut_{ab}) with lenfn going up by exactly 1."""
    out = []
    lu = lenfn(u)
    for a in range(1, n+1):
        for b in range(a+1, n+1):
            v = rmul_t(u, a, b)
            if lenfn(v) == lu+1:
                out.append((a, b, v))
    return out

def chains(u, w, n, lenfn=length):
    """all saturated chains u -> w as lists of (a,b) labels."""
    m = lenfn(w)-lenfn(u)
    if m < 0:
        return
    if u == w:
        yield []
        return
    for (a, b, v) in covers(u, n, lenfn):
        if lenfn(w)-lenfn(v) < 0:
            continue
        for rest in chains(v, w, n, lenfn):
            yield [(a, b)] + rest

def T_tensor(u, w, n, lenfn=length):
    """T_{w/u} as a Counter on label tuples (p_1,...,p_m), p_i in {1..n-1}.
    Coefficient of alpha_{p_1} (x) ... (x) alpha_{p_m}."""
    T = Counter()
    for ch in chains(u, w, n, lenfn):
        # expand each t_{ab} = alpha_a + ... + alpha_{b-1}
        parts = [list(range(a, b)) for (a, b) in ch]
        acc = [()]
        for P in parts:
            acc = [t+(p,) for t in acc for p in P]
        for t in acc:
            T[t] += 1
    return T

def T_free(u, w, n, lenfn=length):
    """The SAME chain sum but in the FREE algebra on transpositions:
    no relation imposed.  Counter on tuples of (a,b) pairs."""
    T = Counter()
    for ch in chains(u, w, n, lenfn):
        T[tuple(ch)] += 1
    return T
