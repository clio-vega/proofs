"""Six-vertex partition function with general weights, and the free-fermion locus.

Vertex convention (t, h, b, r) = (top spin, horizontal-in from left, bottom spin,
horizontal-out to right); spin conservation t + h = b + r.

   (0,0,0,0) a1     (1,1,1,1) a2
   (0,1,0,1) b1     (1,0,1,0) b2
   (1,0,0,1) c1     (0,1,1,0) c2

Free-fermion condition:  a1*a2 + b1*b2 = c1*c2.
"""
from collections import defaultdict
import sympy as sp

SIX = {(0,0,0,0):'a1', (1,1,1,1):'a2', (0,1,0,1):'b1',
       (1,0,1,0):'b2', (1,0,0,1):'c1', (0,1,1,0):'c2'}


def row_transfer(top_vals, left_val, right_val, m, w):
    """All (bottom_state, weight) for one row; w maps vertex name -> weight."""
    out = []
    def step(j, h, bottoms, wt):
        if wt == 0:
            return
        if j == m:
            if h == right_val:
                out.append((tuple(bottoms), wt))
            return
        t = top_vals[j]
        for b in (0, 1):
            for r in (0, 1):
                cfg = (t, h, b, r)
                if cfg in SIX:
                    step(j+1, r, bottoms + [b], wt * w[SIX[cfg]])
    step(0, left_val, [], 1)
    return out


def Z(n, m, top, bot, left, right, weights):
    """weights[i] = dict for row i."""
    states = {tuple(top): sp.Integer(1)}
    for i in range(n):
        new = defaultdict(lambda: sp.Integer(0))
        for st, wv in states.items():
            for ns, rw in row_transfer(st, left[i], right[i], m, weights[i]):
                new[ns] += wv * rw
        states = {k: sp.expand(v) for k, v in new.items()}
    return sp.expand(states.get(tuple(bot), sp.Integer(0)))


def schur_weights(x):
    """The five-vertex free-fermion (Schur) point for a row with parameter x.
    a1=a2=c1=c2=1, b1=x, b2=0.  Free-fermion: 1*1 + x*0 = 1 = 1*1.  OK."""
    return {'a1':sp.Integer(1),'a2':sp.Integer(1),'b1':x,'b2':sp.Integer(0),
            'c1':sp.Integer(1),'c2':sp.Integer(1)}


def ff_weights(x, eps):
    """One-parameter free-fermion deformation off the Schur point.
    Keep a1=a2=c1=1, b1=x, b2=eps; free-fermion forces c2 = a1*a2 + b1*b2 = 1 + x*eps."""
    return {'a1':sp.Integer(1),'a2':sp.Integer(1),'b1':x,'b2':eps,
            'c1':sp.Integer(1),'c2':sp.Integer(1)+x*eps}


def check_ff(w):
    return sp.simplify(w['a1']*w['a2'] + w['b1']*w['b2'] - w['c1']*w['c2']) == 0


if __name__ == "__main__":
    x = sp.symbols('x1:6')
    eps = sp.Symbol('eps')
    print("free-fermion at Schur point :", check_ff(schur_weights(x[0])))
    print("free-fermion after deform   :", check_ff(ff_weights(x[0], eps)))

    print("\n=== locating the Schur boundary condition (control) ===")
    n, m = 2, 5
    ws = [schur_weights(x[i]) for i in range(n)]
    import itertools
    for topset in itertools.combinations(range(m), n):
        top = [1 if j in topset else 0 for j in range(m)]
        for botset in itertools.combinations(range(m), n):
            bot = [1 if j in botset else 0 for j in range(m)]
            z = Z(n, m, top, bot, [0]*n, [0]*n, ws)
            if z != 0:
                print("  top=%s bot=%s  Z = %s" % (topset, botset, sp.factor(z)))
