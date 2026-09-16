"""order.py -- measure ord(N_e) directly.

ord(X) <= d   iff   C^{(d+1)}_{a_1..a_{d+1}} = 0 for all a_i, where

    C^{(d)}_{a_1..a_d} = [...[[X, M_{p_a1}], M_{p_a2}], ..., M_{p_ad}]
                       = sum_{S subset [d]} (-1)^{|S|} M_{p_{a_S}} X M_{p_{a_{S^c}}}

(the p's commute, so only the SET matters).

Harness rules obeyed:
 (1) sweep the a-tuples with itertools.product, NOT combinations -- the tuple
     (a_1,...,a_d) is ordered but the operator is symmetric in the a_i; we
     nevertheless sweep the full product so the quantifier in the write-up is
     "for all ordered tuples in {1..A}^d".
 (2) report C^{(d)} nonzero on the SAME data on which C^{(d+1)} vanishes.
 (3) vary the source state lam over many shapes, including the vacuum.
"""
import itertools, sys
from fast import op_p, op_N, op_Psi

def parts(n, cap=None):
    if cap is None: cap = n
    if n == 0: return [()]
    return [(p,) + r for p in range(min(n, cap), 0, -1) for r in parts(n - p, p)]

def mult(state, aa):
    for a in aa: state = op_p(state, a)
    return state

def nested(X, e, aa, lam):
    """C^{(d)}_{aa} applied to |lam>."""
    d = len(aa); out = {}
    for r in range(d + 1):
        for S in itertools.combinations(range(d), r):
            sgn = (-1)**r
            Sc = [aa[i] for i in range(d) if i not in S]
            st = mult({lam: 1}, Sc)
            st = X(st, e)
            st = mult(st, [aa[i] for i in S])
            for k, v in st.items():
                out[k] = out.get(k, 0) + sgn * v
    return {k: v for k, v in out.items() if v}

def scan(X, name, es, A, ds, lams, verbose=True):
    for e in es:
        row = []
        for d in ds:
            nz = 0; tot = 0; witness = None
            for aa in itertools.product(range(1, A + 1), repeat=d):
                for lam in lams:
                    tot += 1
                    r = nested(X, e, aa, lam)
                    if r:
                        nz += 1
                        if witness is None: witness = (aa, lam, dict(list(r.items())[:3]))
            row.append((d, nz, tot, witness))
        print('%s  e=%d' % (name, e))
        for (d, nz, tot, w) in row:
            print('    C^(%d): %5d/%5d nonzero' % (d, nz, tot),
                  ('  witness a=%s lam=%s -> %s' % w) if w else '   ==> VANISHES')
        sys.stdout.flush()

if __name__ == '__main__':
    lams = [()] + parts(1) + parts(2) + parts(3) + parts(4) + parts(5) + parts(6)
    scan(op_N, 'N_e', [2, 3, 4], 3, [1, 2, 3], lams)

def deep(X, e, dmax, lams, A=1):
    for d in range(1, dmax + 1):
        nz = 0; tot = 0; w = None
        for aa in itertools.product(range(1, A + 1), repeat=d):
            for lam in lams:
                tot += 1
                r = nested(X, e, aa, lam)
                if r:
                    nz += 1
                    if w is None: w = (aa, lam, dict(sorted(r.items())[:2]))
        print('  e=%d  C^(%d): %d/%d nonzero' % (e, d, nz, tot),
              ('  witness %s %s -> %s' % w) if w else '  ==> VANISHES')
        sys.stdout.flush()
