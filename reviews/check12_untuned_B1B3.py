"""check12_untuned_B1B3.py -- an UNTUNED extra detector for Lyra's level-2 leg.

Uglov math/9905196: B_m := sum_i X_i^m lies in Z(affH), the CENTRE of the
affine Hecke algebra, which Bernstein's theorem identifies with symmetric
Laurent polynomials in the commuting X_i.  A centre is commutative, so
    [B_{-a}, B_{-b}] = 0   for ALL a,b > 0,
not merely for the pair (1,2) that Lyra tested.  (Her stated reason -- the
Heisenberg relation "[B_m,B_m']=0 unless m+m'=0" -- is a different algebra;
the correct warrant is commutativity of Z(affH), and it is STRONGER.)

(1,3) and (2,3) were never used to fix any choice in either implementation,
so they are free, untuned tests of both straighteners.
"""
import sys
import sympy as sp

sys.path.insert(0, "/home/clio/projects/probes/2026-08-31-route1-diff")
sys.path.insert(0, "/tmp/lyra-math/level2-route1")
import route3_uglov as MINE
import wedge as LYRA
import heisenberg as LYRA_H

def partitions_upto(nmax):
    out = [()]
    def rec(rem, cap, cur):
        for part in range(min(rem, cap), 0, -1):
            nc = cur + (part,)
            out.append(nc)
            if rem - part > 0:
                rec(rem - part, part, nc)
    for t in range(1, nmax + 1):
        rec(t, t, ())
    seen, uniq = set(), []
    for p in out:
        if p not in seen:
            seen.add(p); uniq.append(p)
    return uniq

def wedge_of(lam, R, s=0):
    p = list(lam) + [0] * (R - len(lam))
    return tuple(p[r] - r + s for r in range(R))

def my_B(state, m, n, ell):
    shift = m * n * ell
    sh = {}
    for idx, c in state.items():
        for r in range(len(idx)):
            new = list(idx); new[r] += shift
            sh[tuple(new)] = sh.get(tuple(new), sp.Integer(0)) + c
    return MINE.straighten(sh, n, ell)

def comm(state, a, b, n, ell, impl):
    x = impl(impl(state, b, n, ell), a, n, ell)
    y = impl(impl(state, a, n, ell), b, n, ell)
    out = {}
    for k in set(x) | set(y):
        d = sp.cancel(sp.together(x.get(k, 0) - y.get(k, 0)))
        if sp.simplify(d) != 0:
            out[k] = d
    return out

if __name__ == "__main__":
    ELL, NMAX, RX, CH, E = 2, 4, 4, [0, 1], [2, 3]
    lams = partitions_upto(NMAX)
    print("=" * 72)
    print("UNTUNED DETECTORS  [B_-1,B_-3]=0  and  [B_-2,B_-3]=0   at ell=2")
    print(f"  e in {E}, |lam|<={NMAX}, R=len(lam)+{RX}, charges {CH}")
    print("=" * 72)
    for (a, b) in [(1, 3), (2, 3)]:
        for label, impl in [("CLIO  route3_uglov", my_B),
                            ("LYRA  level2-route1", LYRA_H.apply_Bminus)]:
            nz = nf = 0; ex = []
            for lam in lams:
                R = len(lam) + RX
                for s in CH:
                    for e in E:
                        st = {wedge_of(lam, R, s): sp.Integer(1)}
                        c = comm(st, a, b, e, ELL, impl)
                        if c:
                            nf += 1
                            if len(ex) < 2:
                                ex.append((e, lam, s, list(c.items())[:2]))
                        else:
                            nz += 1
            print(f"  [B_-{a},B_-{b}]  {label}:  ZERO {nz}/{nz+nf}"
                  + (f"   FAILURES {nf}" if nf else ""))
            for z in ex:
                print("      ", z)
    print("=" * 72)
