"""C4 -- REFUSAL CONTROL ON THE CITATION, not on the statement (Q).

The one step of the chain I cannot re-read this session is Branden-Huh
\label{CorollaryConvolution} (l.2416), recorded as

      N(f), N(g) Lorentzian  =>  N(fg) Lorentzian.

The recorded form does not say what the hypotheses are beyond that.  So: try to
REFUTE the recorded form in 3 variables.  If a counterexample exists the recorded
form is mis-transcribed and the chain is dead; if none exists over an exhaustive
small range plus a large randomised range, the recorded form is at least not
mis-transcribed in a way that matters at 3 variables.

Also recorded for each surviving pair: whether nu_f = log c is M-concave, so that
I can see whether the corollary is being exercised BEYOND the reach of
\label{normalizedcoefficients} (if every Lorentzian N(f) in range had nu_f
M-concave, the test would be vacuous as a test of the convolution corollary).
"""
import random, itertools, sys
from lor import *

def mono(d, n=3):
    def rec(k, rem):
        if k == 1: yield (rem,); return
        for v in range(rem+1):
            for t in rec(k-1, rem-v): yield (v,)+t
    return list(rec(n, d))

def is_M_concave(c):
    """nu(alpha)=log c_alpha M-concave:  dom M-convex and for all alpha,beta in dom,
       all i with alpha_i>beta_i, exists j with alpha_j<beta_j s.t. both exchanged
       points lie in dom and  c_alpha c_beta <= c_{a-ei+ej} c_{b+ei-ej}."""
    supp = [a for a,v in c.items() if v != 0]
    S = set(supp)
    if not is_Mconvex(supp): return False
    n = len(supp[0])
    for a in supp:
        for b in supp:
            for i in range(n):
                if a[i] > b[i]:
                    ok = False
                    for j in range(n):
                        if a[j] < b[j]:
                            a2=list(a); a2[i]-=1; a2[j]+=1; a2=tuple(a2)
                            b2=list(b); b2[i]+=1; b2[j]-=1; b2=tuple(b2)
                            if a2 in S and b2 in S and c[a]*c[b] <= c[a2]*c[b2]:
                                ok=True; break
                    if not ok: return False
    return True

# ---- step 1: harvest all f of degree 2 and 3 in 3 vars with small coeffs and N(f) Lorentzian
def harvest(d, maxc, cap=None, randomised=0):
    ms = mono(d)
    out = []
    if randomised:
        for _ in range(randomised):
            c = {a: random.randint(0, maxc) for a in ms}
            c = {a:v for a,v in c.items() if v}
            if not c: continue
            if is_N_lorentzian(c): out.append(c)
    else:
        for vals in itertools.product(range(maxc+1), repeat=len(ms)):
            c = {a:v for a,v in zip(ms,vals) if v}
            if not c: continue
            if is_N_lorentzian(c): out.append(c)
            if cap and len(out) >= cap: break
    return out

random.seed(1004)
print("harvesting Lorentzian N(f) in 3 variables ...")
H2 = harvest(2, 3)                       # exhaustive, 4^6 = 4096 candidates
H3 = harvest(3, 2)                       # exhaustive, 3^10 = 59049 candidates
print(f"  degree 2, coeffs 0..3 : {len(H2)} Lorentzian out of {4**6}")
print(f"  degree 3, coeffs 0..2 : {len(H3)} Lorentzian out of {3**10}")
nmc2 = sum(1 for c in H2 if not is_M_concave(c))
nmc3 = sum(1 for c in H3 if not is_M_concave(c))
print(f"  of these, nu_f NOT M-concave: degree 2: {nmc2}/{len(H2)}, degree 3: {nmc3}/{len(H3)}")
print("  (so the test is not vacuous: it exercises pairs outside normalizedcoefficients)")

def test_pairs(A, B, label, limit=None):
    tot = fail = beyond = 0
    fails = []
    pairs = [(f,g) for f in A for g in B]
    if limit and len(pairs) > limit:
        pairs = random.sample(pairs, limit)
    for f,g in pairs:
        tot += 1
        if not (is_M_concave(f) and is_M_concave(g)): beyond += 1
        h = poly_mult(f,g)
        ok, w = is_N_lorentzian(h, return_witness=True)
        if not ok:
            fail += 1
            if len(fails) < 3: fails.append((f,g,w))
    print(f"  {label}: pairs enumerated {tot}, counterexamples {fail}, "
          f"pairs where at least one factor is beyond normalizedcoefficients: {beyond}")
    for f,g,w in fails: print("     COUNTEREXAMPLE", f, g, w)
    return fail

print("\ntesting the recorded form  N(f),N(g) Lorentzian => N(fg) Lorentzian:")
tf = 0
tf += test_pairs(H2, H2, "deg 2 x deg 2 -> deg 4 (exhaustive)")
tf += test_pairs(H2, H3, "deg 2 x deg 3 -> deg 5", limit=4000)
tf += test_pairs(H3, H3, "deg 3 x deg 3 -> deg 6", limit=3000)
print("\nTOTAL counterexamples to the recorded form:", tf)
