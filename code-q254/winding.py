"""Q254: does the tested range actually WIND?  (Q253 lesson: calibrate in the
degenerate regime, then make sure the tested range leaves it.)

A cylindric shape in the bead model is a tuple x_1<...<x_m<x_1+n.  The ONLY
place the cylinder enters the tableau recursion is the wrap constraint
   rho_m < nu_1 + n      (cyl.hstrip, i=m-1, "nxt = nu[0]+n").
If n is replaced by a large N (keeping the same tuple), that constraint stops
binding and the model computes the ORDINARY skew Schur function of the
unrolled shape.  So:

   the shape WINDS at ell variables  <=>  W(n,m,mu,lam,ell) != W(N,m,mu,lam,ell)
                                          for N >> 0.

Also computed here: the monomial multiplicities and the Schur expansion, to
test whether s^c is Schur-positive.
"""
import sys, itertools
from collections import defaultdict
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
import cyl as C   # vendored copy, see cyl.py header
from rigidity import dominates, sort_part, partitions_le


def weights_counted(mu, lam, n, m, ell):
    """monomial multiplicities: weight tuple -> #cylindric SSYT"""
    S = C.interval_shapes(mu, lam, n, m)
    idx = {s: j for j, s in enumerate(S)}
    adj = [[] for _ in S]
    for a, sa in enumerate(S):
        for b, sb in enumerate(S):
            if C.hstrip(sa, sb, n, m):
                adj[a].append((b, C.size(sa, sb)))
    states = {idx[mu]: {(): 1}}
    for t in range(ell):
        nxt = {}
        for a, ws in states.items():
            for b, sz in adj[a]:
                d = nxt.setdefault(b, defaultdict(int))
                for w, c in ws.items():
                    d[w + (sz,)] += c
        states = nxt
    return dict(states.get(idx[lam], {}))


def winds(mu, lam, n, m, ell, big=None):
    if big is None:
        big = max(n, lam[-1] - mu[0] + 2) + n + 4
    W1 = set(C.weights(mu, lam, n, m, ell))
    W2 = set(C.weights(mu, lam, big, m, ell))
    return W1 != W2, W1, W2


# ---- monomial -> Schur, via Kostka inverse on partitions of N with <=ell parts
def kostka(lam, mu):
    """K_{lam,mu} = #SSYT of shape lam and content mu, by brute force."""
    lam = list(lam)
    cells = [(i, j) for i in range(len(lam)) for j in range(lam[i])]
    cnt = 0
    mu = list(mu)
    def rec(k, T, rem):
        nonlocal cnt
        if k == len(cells):
            cnt += 1; return
        i, j = cells[k]
        lo = 1
        if j > 0: lo = max(lo, T[(i, j-1)])
        if i > 0: lo = max(lo, T[(i-1, j)] + 1)
        for v in range(lo, len(mu)+1):
            if rem[v-1] == 0: continue
            T[(i, j)] = v; rem[v-1] -= 1
            rec(k+1, T, rem)
            rem[v-1] += 1; del T[(i, j)]
    rec(0, {}, mu[:])
    return cnt


def schur_expansion(mcoeffs, N, ell):
    """mcoeffs: weight tuple -> multiplicity.  Returns dict partition -> coeff."""
    parts = [p for p in partitions_le(N, ell)]
    # monomial coefficients m_lam
    mono = {}
    for w, c in mcoeffs.items():
        p = sort_part(w)
        mono[p] = c            # all rearrangements have the same coefficient
    # f = sum_p mono[p] m_p ;  s_lam = sum_mu K_{lam,mu} m_mu  => solve upper-triangular
    parts_sorted = sorted(parts, key=lambda p: (-sum(1 for _ in p), p), reverse=True)
    # order by dominance: process dominance-maximal first
    order = sorted(parts, key=lambda p: tuple(itertools.accumulate(list(p)+[0]*(ell-len(p)))),
                   reverse=True)
    rem = dict(mono)
    out = {}
    for lam in order:
        c = rem.get(lam, 0)
        if c == 0: continue
        out[lam] = c
        for mu in parts:
            k = kostka(lam, mu)
            if k:
                rem[mu] = rem.get(mu, 0) - c * k
    assert all(v == 0 for v in rem.values()), rem
    return out


if __name__ == "__main__":
    from cylindric import cyl_shapes
    tot = wind = 0
    nonpos = []
    examples = []
    for n in range(2, 6):
        for m in range(1, n):
            for mu in cyl_shapes(n, m):
                for lam in itertools.product(*[range(mu[i], mu[i]+n+1) for i in range(m)]):
                    if not C.is_shape(lam, n, m) or not C.contains(lam, mu): continue
                    N = C.size(mu, lam)
                    if N == 0 or N > 6: continue
                    for ell in range(1, 5):
                        W = set(C.weights(mu, lam, n, m, ell))
                        if not W: continue
                        tot += 1
                        w, W1, W2 = winds(mu, lam, n, m, ell)
                        if w:
                            wind += 1
                            if len(examples) < 3: examples.append((n,m,mu,lam,ell,sorted(W1),sorted(W2)))
                        if w and N <= 5 and ell <= 4:
                            mc = weights_counted(mu, lam, n, m, ell)
                            try:
                                sx = schur_expansion(mc, N, ell)
                            except AssertionError:
                                continue
                            if any(v < 0 for v in sx.values()):
                                nonpos.append((n, m, mu, lam, ell, sx))
    print(f"WINDING: {wind}/{tot} (shape,ell) instances have a support that CHANGES when the")
    print(f"         wrap constraint is relaxed -- i.e. the cylinder is doing work.")
    for e in examples:
        print(f"   e.g. n={e[0]} m={e[1]} mu={e[2]} lam={e[3]} ell={e[4]}")
        print(f"        cylindric supp {e[5]}")
        print(f"        unrolled  supp {e[6]}")
    print()
    print(f"SCHUR-POSITIVITY: {len(nonpos)} winding instances with a NEGATIVE Schur coefficient")
    for e in nonpos[:4]:
        print(f"   n={e[0]} m={e[1]} mu={e[2]} lam={e[3]} ell={e[4]}: {e[5]}")
