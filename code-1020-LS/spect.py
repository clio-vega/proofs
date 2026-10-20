r"""The GENERAL form: ANY set of spectators disjoint from the row index.

CLAIM: let Sigma subset [N] with Sigma cap A = empty.  Pin w_s = c_s, arbitrary
finite nonzero constants, for s in Sigma; put w_c = T^{-e_c} (e strictly
increasing on [N]\Sigma) for the rest; q = 0.  Then

      Atilde^X_{A,B} ---> [A >= B]   as T -> 0+,   X in {B,C}.

Rick's Lemma L is the case |Sigma| = 2 with c_a = x, c_b = 1/x generic;
hyp:L^S is the case |Sigma| = 1 with c_m = eps in {+1,-1}.  Same lemma.

DOMAIN: N = 1..4; X in {B,C}; EVERY Sigma subset [N] (including empty and
        non-singletons); every (A,B) with A cap Sigma = empty.
        Spectator values drawn from a deliberately ADVERSARIAL list including
        +1, -1, and Rick's reciprocal pair (x, 1/x) -- NOT generic.
        T = 1/10, 1/100, 1/1000, exact rationals, decay-exponent test.
"""
import sys, math, itertools
import sympy as sp
sys.path.insert(0,'.'); import q484_engine as E
E.q = sp.Integer(0)

# adversarial spectator values: the two roots of unity, and a reciprocal pair
SPEC = [sp.Integer(1), sp.Integer(-1), sp.Rational(5,3), sp.Rational(3,5)]
Ts = [sp.Rational(1,10), sp.Rational(1,100), sp.Rational(1,1000)]

def ok_decay(errs):
    if all(x == 0 for x in errs): return True
    ps = []
    for i in range(len(errs)-1):
        if errs[i] == 0 or errs[i+1] == 0: ps.append(float('inf')); continue
        ps.append(math.log(float(errs[i]/errs[i+1]))/math.log(float(Ts[i]/Ts[i+1])))
    return all(p > 0.5 for p in ps)

good = plant = bad = 0; fails = []
for N in range(1,5):
    P = E.subsets(N)
    for Sigma in P:                                    # every spectator set
        nonspec = [c for c in range(1,N+1) if c not in Sigma]
        e = {c: i+1 for i, c in enumerate(nonspec)}    # strictly increasing
        cs = {s: SPEC[i % len(SPEC)] for i, s in enumerate(Sigma)}
        for X in ('C','B'):
            for A in P:
                if set(A) & set(Sigma): continue
                for B in P:
                    pred = sp.Integer(1) if E.geq(A,B) else sp.Integer(0)
                    errs = []
                    for T in Ts:
                        w = [cs[c] if c in Sigma else sp.Rational(1,1)/T**e[c]
                             for c in range(1,N+1)]
                        v = sp.Rational(sp.cancel(E.Atilde(N,w,A,B,X=X,simplify=False)))
                        errs.append(abs(v - pred))
                    if ok_decay(errs):
                        if pred == 1: good += 1
                        else: plant += 1
                    else:
                        bad += 1
                        if len(fails)<8: fails.append((N,X,Sigma,A,B,int(pred),[str(x) for x in errs]))
print(f"N<=4, every spectator set Sigma, both types, adversarial spectator values:")
print(f"  A>=B  (limit 1): {good}")
print(f"  A not>=B (plant, limit 0): {plant}")
print(f"  TOTAL {good+plant+bad}   FAILURES {bad}")
for f in fails: print("   ", f)
