"""First-hand verification of Rick's Remarks A and B (peer-claimed -> mine only if I derive them).

Remark A.  Closed form  g^t_i = min(lambda_i, mu_{i+t} - t), hence
           l_0 = min{t : mu_{i+t} - t >= lambda_i for all i} < infinity for EVERY mu subset lambda,
           so T_l != empty  iff  l >= l_0  unconditionally, and the l_min = infinity
           branch of note 3 Lemma 4.1 is VACUOUS.

  My derivation (to be tested, not assumed):
    unrolling g^t_i = min(lambda_i, g^{t-1}_{i+1} - 1) gives
        g^t_i = min( min_{0<=s<t} (lambda_{i+s} - s),  mu_{i+t} - t ).
    Rick's two-term form needs the intermediate terms ABSORBED:
        lambda_{i+s} - s >= lambda_i   for all s >= 0.                        (*)
    (*) holds because bead positions are STRICTLY increasing: lambda_{i+1} >= lambda_i + 1.
    This is the direction the brief flagged -- confirm it, do not assume it.
    Finiteness: mu_{i+km} - km = mu_i + k(n-m) -> infinity since m < n, so t = km works
    for k large.  No induction on lambda needed for THIS part.

Remark B.  gamma = (|g^t| - |g^{t-1}|)_t is already weakly decreasing, so gamma = lambdahat.
"""
import sys, itertools
from collections import Counter
sys.path.insert(0, '/home/clio/projects/peers/rick/proofs')
import types
src = open('/home/clio/projects/peers/rick/proofs/2026-09-25-clio288-greedy-check.py').read().split('stats = Counter()')[0]
rick = types.ModuleType('rick'); exec(compile(src, 'rick', 'exec'), rick.__dict__)
ext, greedy, shapes, srt = rick.ext, rick.greedy, rick.shapes, rick.srt

def main(NMAX=7, DMAX=9):
    out = Counter(); worst_k = 0
    for n in range(2, NMAX+1):
        for m in range(1, n):
            for mu in shapes(n, m, 0, n):
                if mu[0] != 0: continue
                for lam in shapes(n, m, 0, n+DMAX+1):
                    if not all(l >= u for l, u in zip(lam, mu)): continue
                    if sum(lam)-sum(mu) > DMAX: continue
                    out['shapes'] += 1
                    # (*) strictly increasing beads => lambda_{i+s} - s >= lambda_i
                    for i in range(m):
                        for s in range(0, 3*m+4):
                            if ext(lam, n, i+1+s) - s < lam[i]:
                                out['absorb_FAIL'] += 1
                                if out['absorb_FAIL'] <= 3: print("  (*) FAIL", n, lam, i, s)
                    g = greedy(mu, lam, n, 400)
                    # closed form, two-term (Rick) and full-min (mine): both vs the recursion
                    for t in range(1, 40):
                        for i in range(m):
                            rec  = g[t][i]
                            rick_form = min(lam[i], ext(mu, n, i+1+t) - t)
                            full = min(min(ext(lam, n, i+1+s) - s for s in range(t)),
                                       ext(mu, n, i+1+t) - t)
                            if rec != full: out['fullmin_FAIL'] += 1
                            if rec != rick_form:
                                out['closedform_FAIL'] += 1
                                if out['closedform_FAIL'] <= 3:
                                    print("  closed form FAIL", n, mu, lam, t, i, rec, rick_form)
                    # finiteness of l_0, and the explicit bound t = k*m
                    l0 = next((t for t, x in enumerate(g) if x == lam), None)
                    if l0 is None:
                        out['l0_INFINITE'] += 1; print("  l0 INFINITE", n, mu, lam); continue
                    out['l0_finite'] += 1
                    k = -(-l0 // m) if m else 0
                    worst_k = max(worst_k, k)
                    # the a-priori bound: t = k*m with mu_i + k(n-m) >= lam_i for all i
                    kb = max((-(-(lam[i]-mu[i]) // (n-m)) for i in range(m)), default=0)
                    if l0 > kb*m:
                        out['bound_FAIL'] += 1
                        if out['bound_FAIL'] <= 3: print("  bound FAIL", n, mu, lam, l0, kb*m)
                    # T_l nonempty iff l >= l0 is Rick's script's own check; here: gamma sorted
                    gam = tuple(sum(g[t])-sum(g[t-1]) for t in range(1, max(l0,1)+1))
                    if any(gam[t] < gam[t+1] for t in range(len(gam)-1)):
                        out['gamma_NOT_SORTED'] += 1
                        if out['gamma_NOT_SORTED'] <= 3: print("  gamma not sorted", n, mu, lam, gam)
                    elif srt(gam) != tuple(x for x in gam if x > 0):
                        out['gamma_ne_lambdahat'] += 1
    print(dict(out)); print("worst ceil(l0/m) =", worst_k)
    bad = sum(v for k, v in out.items() if 'FAIL' in k or 'INFINITE' in k or 'NOT_SORTED' in k or '_ne_' in k)
    print("REMARKS A,B ALL PASS" if bad == 0 else f"{bad} FAILURES")
    return bad == 0

if __name__ == "__main__":
    a=[int(x) for x in sys.argv[1:]] or [7,9]
    sys.exit(0 if main(*a) else 1)
