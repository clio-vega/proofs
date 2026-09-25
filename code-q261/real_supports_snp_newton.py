"""
V4: T1 (SNP) and T2 (Newton = P_lambdahat) tested on the REAL cylindric
supports W_ell, obtained from Rick's independent enumerator
peers/rick/proofs/2026-09-25-clio288-greedy-check.py, which builds chains from
the RAW strip condition (1) and never touches Lemma 3.1 or my greedy chain.

Deliberately NOT routed through thm:main.  The bounding inequalities are read
off W itself:
      M_r := max_{alpha in W} (sum of the r largest entries of alpha),
so  conv(W) is contained in Q_W = {x >= 0 : x([l]) = d, x(S) <= M_{|S|}}
with no reference to lambdahat.  T1 is then the assertion Q_W cap Z^l = W.

T2 is certified CONSTRUCTIVELY: for each alpha in W we run the Robin Hood
chain of Lemma 4.1 from lambdahat down to sort(alpha) and produce an exact
rational convex combination of PERMUTATIONS of lambdahat equal to alpha.
No LP, no floating point -- an explicit certificate.
"""
import sys, itertools
from fractions import Fraction
from collections import Counter
sys.path.insert(0, '/home/clio/projects/peers/rick/proofs')
sys.path.insert(0, '/home/clio/projects/proofs/code-q261')

# import Rick's primitives without running his main loop
import types
src = open('/home/clio/projects/peers/rick/proofs/2026-09-25-clio288-greedy-check.py').read()
src = src.split('stats = Counter()')[0]
rick = types.ModuleType('rick'); rick.__dict__['__name__'] = 'rick'
exec(compile(src, 'rick', 'exec'), rick.__dict__)

from rado_from_robin_hood import sort_desc, psums, dominates, compositions, unroll, chain

def topr(alpha, r):
    return sum(sorted(alpha, reverse=True)[:r])

def check_instance(mu, lam, n, ell, W, lamhat, out):
    d = sum(lam) - sum(mu)
    Wset = set(W)
    # ---- inequalities read off W itself, NOT from lambdahat
    M = [0] + [max(topr(a, r) for a in Wset) for r in range(1, ell + 1)]
    # sanity: W really does satisfy them, and they are the tightest such
    for a in Wset:
        assert sum(a) == d
        for r in range(1, ell + 1):
            assert topr(a, r) <= M[r]
    # ---- T1: Q_W cap Z^l = W
    QW = set()
    for a in compositions(d, ell):
        if all(sum(a[i] for i in S) <= M[len(S)]
               for k in range(ell + 1) for S in itertools.combinations(range(ell), k)):
            QW.add(a)
    if QW != Wset:
        out['T1_fail'] += 1
        if out['T1_fail'] <= 3:
            print("  T1 FAIL", n, mu, lam, ell, "extra:", sorted(QW - Wset)[:4],
                  "missing:", sorted(Wset - QW)[:4])
        return
    out['T1_ok'] += 1
    out['T1_points'] += len(QW)
    # ---- independent: M_r should equal Lambda_r of lambdahat  (a consequence, checked not assumed)
    lhpad = tuple(lamhat) + (0,) * (ell - len(lamhat))
    Lam = [0] + psums(lhpad)
    if M[:ell + 1] != Lam[:ell + 1]:
        out['M_ne_Lambda'] += 1
        print("  M != Lambda", n, mu, lam, ell, M, Lam)
    # ---- T2: constructive certificate for every alpha in W
    for a in Wset:
        sigma = tuple(sorted(a, reverse=True))
        lh = tuple(lamhat) + (0,) * (ell - len(lamhat))
        if not dominates(lh, sigma):
            out['T2_fail'] += 1; print("  T2 dominance FAIL", lamhat, sigma); continue
        chain(lh, sigma)                       # obligations: gap >= 2, Phi decreasing
        combo = unroll(lh, sigma)
        if (any(c < 0 for c in combo.values()) or sum(combo.values()) != 1
                or any(sort_desc(p) != sort_desc(lh) for p in combo)):
            out['T2_fail'] += 1; print("  T2 certificate FAIL", lamhat, sigma); continue
        got = tuple(sum(c * p[k] for p, c in combo.items()) for k in range(ell))
        if got != tuple(Fraction(x) for x in sigma):
            out['T2_fail'] += 1; print("  T2 wrong point", lamhat, sigma, got); continue
        out['T2_certs'] += 1
    # ---- P_lambdahat subset conv(W): every permutation of lambdahat is IN W
    for p in set(itertools.permutations(tuple(lamhat) + (0,) * (ell - len(lamhat)))):
        if p not in Wset:
            out['perm_fail'] += 1
            if out['perm_fail'] <= 3: print("  PERM FAIL", n, mu, lam, ell, p)

def main(NMAX=5, DMAX=6, LMAX=4):
    out = Counter()
    for n in range(2, NMAX + 1):
        for m in range(1, n):
            for mu in rick.shapes(n, m, 0, n):
                if mu[0] != 0: continue
                for lam in rick.shapes(n, m, 0, n + DMAX + 1):
                    if not all(l >= u for l, u in zip(lam, mu)): continue
                    d = sum(lam) - sum(mu)
                    if d > DMAX: continue
                    g = rick.greedy(mu, lam, n, 200)
                    l0 = next((t for t, x in enumerate(g) if x == lam), None)
                    if l0 is None: out['l0_inf'] += 1; continue
                    for ell in range(1, LMAX + 1):
                        W = rick.weights(mu, lam, n, ell)
                        if not W: continue
                        out['instances'] += 1
                        gam = tuple(sum(g[t]) - sum(g[t - 1]) for t in range(1, ell + 1))
                        lamhat = rick.srt(gam)
                        check_instance(mu, lam, n, ell, W, lamhat, out)
    print(dict(out))
    ok = not (out['T1_fail'] or out['T2_fail'] or out['perm_fail'] or out['M_ne_Lambda'])
    print("V4 ALL PASS" if ok else "V4 FAILURES")
    return ok

if __name__ == "__main__":
    a = [int(x) for x in sys.argv[1:]] or [5, 6, 4]
    sys.exit(0 if main(*a) else 1)
