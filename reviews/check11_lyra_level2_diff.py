"""check11_lyra_level2_diff.py -- Clio's INDEPENDENT Route-2 crosscheck of
Lyra's level2-route1 at ell=2  (peer review 2026-09-02).

Answers Lyra's ask (UID 683): "does your Route-2 straightener diff GREEN
against level2-route1 at ell=2?"

My side uses MY OWN straightener, probes/2026-08-31-route1-diff/route3_uglov.py,
unchanged, plus a Heisenberg wrapper I wrote from the PRIMARY source --
Uglov math/9905196 eq. (e:Bo):
    B_m(u_k) = sum_j u_{k_1} ^ ... ^ u_{k_j - n*l*m} ^ ... ^ u_{k_r}
so B_{-m} shifts each bead by +m*n*l = +m*N.  This is Uglov's own definition,
NOT a value tuned to make a detector pass.

Comparison layer uses sp.cancel (my straightener uses sp.expand internally,
which does NOT cancel the R4 /(q+q^{-1}) denominators -- Lyra's bug #1).
"""
import sys, itertools
import sympy as sp

sys.path.insert(0, "/home/clio/projects/probes/2026-08-31-route1-diff")
sys.path.insert(0, "/tmp/lyra-math/level2-route1")

import route3_uglov as MINE          # Clio, Route 2
import wedge as LYRA                 # Lyra, level2-route1
import heisenberg as LYRA_H

q = sp.Symbol("q")


def partitions_upto(nmax):
    out = [()]
    def rec(remaining, cap, cur):
        for part in range(min(remaining, cap), 0, -1):
            nc = cur + (part,)
            out.append(nc)
            if remaining - part > 0:
                rec(remaining - part, part, nc)
    for total in range(1, nmax + 1):
        rec(total, total, ())
    seen, uniq = set(), []
    for p in out:
        if p not in seen:
            seen.add(p); uniq.append(p)
    return uniq


def wedge_of(lam, R, s=0):
    padded = list(lam) + [0] * (R - len(lam))
    return tuple(padded[r] - r + s for r in range(R))


def my_Bminus(state, m, n, ell):
    """Uglov (e:Bo): B_{-m} shifts each bead by +m*N, N = n*ell."""
    shift = m * n * ell
    shifted = {}
    for idx, coeff in state.items():
        for r in range(len(idx)):
            new = list(idx); new[r] = new[r] + shift
            key = tuple(new)
            shifted[key] = shifted.get(key, sp.Integer(0)) + coeff
    return MINE.straighten(shifted, n, ell)


def norm(d):
    """Canonicalise a {tuple->coeff} dict with cancel; drop true zeros."""
    out = {}
    for k, v in d.items():
        c = sp.cancel(sp.together(v))
        if sp.simplify(c) != 0:
            out[k] = sp.cancel(c)
    return out


def diff(a, b):
    a, b = norm(a), norm(b)
    bad = {}
    for k in set(a) | set(b):
        d = sp.simplify(sp.cancel(a.get(k, 0) - b.get(k, 0)))
        if d != 0:
            bad[k] = d
    return bad


if __name__ == "__main__":
    E, ELL, NMAX, RX, CH = [2, 3], 2, 5, 4, [0, 1]
    lams = partitions_upto(NMAX)
    print("=" * 72)
    print("CLIO Route-2  vs  LYRA level2-route1 @5ab2c33   ell=2  B_{-1} DIFF")
    print(f"  e in {E}, ell={ELL}, |lam|<={NMAX}, R=len(lam)+{RX}, charges {CH}")
    print(f"  configurations = {len(lams)*len(CH)*len(E)}")
    print("=" * 72)

    ncfg = nagree = ndis = 0
    disagreements = []
    for lam in lams:
        R = len(lam) + RX
        for s in CH:
            u = wedge_of(lam, R, s)
            for e in E:
                ncfg += 1
                st = {u: sp.Integer(1)}
                mine = my_Bminus(st, 1, e, ELL)
                hers = LYRA_H.apply_Bminus(st, 1, e, ELL)
                d = diff(mine, hers)
                if d:
                    ndis += 1
                    if len(disagreements) < 5:
                        disagreements.append((e, lam, s, d))
                else:
                    nagree += 1
    print(f"\n[B_-1 DIFF]  agree {nagree}/{ncfg}   disagree {ndis}")
    for (e, lam, s, d) in disagreements:
        items = list(d.items())[:3]
        print(f"   e={e} lam={lam} s={s}: " +
              "; ".join(f"{k}:{v}" for k, v in items))

    # my own independent detector, using MY straightener only
    print("\n[CLIO-SIDE DETECTOR]  [B_-1,B_-2]=0 using MY straightener:")
    nz = nf = nt = 0
    fails = []
    for lam in lams:
        R = len(lam) + RX
        for s in CH:
            u = wedge_of(lam, R, s)
            for e in E:
                nt += 1
                st = {u: sp.Integer(1)}
                b1b2 = my_Bminus(my_Bminus(st, 2, e, ELL), 1, e, ELL)
                b2b1 = my_Bminus(my_Bminus(st, 1, e, ELL), 2, e, ELL)
                comm = {}
                for k in set(b1b2) | set(b2b1):
                    comm[k] = b1b2.get(k, 0) - b2b1.get(k, 0)
                comm = norm(comm)
                if comm:
                    nf += 1
                    if len(fails) < 3:
                        fails.append((e, lam, s, list(comm.items())[:2]))
                else:
                    nz += 1
    print(f"    ZERO on {nz}/{nt}  (nonzero {nf})")
    for f in fails:
        print("     ", f)
    print("\n" + "=" * 72)
