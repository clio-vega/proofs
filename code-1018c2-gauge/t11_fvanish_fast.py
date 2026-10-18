"""Lemma (F-vanishing), checked by DIRECT VANISHING at exact rational points of D.

t8 computed ord_D(F) symbolically and did not return in an hour.  The claim is
only ord >= 1, i.e. F = 0 ON the divisor, so evaluate there: on
D = {1 - q w_c/w_e = 0} the function field is C(w) with q = w_e/w_c, so set
q := w_e/w_c at exact rational w and ask whether F vanishes.  Exact arithmetic
in Q; three generic seeds.

Claim:  e < c  =>  F = 0 on D whenever  e in M
        e > c  =>  F = 0 on D whenever  {e,c} subset M
PLANTS (must NOT vanish): the same test when the hypothesis fails.
"""
import sympy as sp, itertools
from q484_engine import q, wvars, F

SEEDS = [[sp.Rational(3,5), sp.Rational(7,11), sp.Rational(13,17)],
         [sp.Rational(4,9), sp.Rational(11,7), sp.Rational(23,19)],
         [sp.Rational(6,5), sp.Rational(17,13), sp.Rational(31,29)]]
def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

score={}
def rec(tag, ok):
    score.setdefault(tag,[0,0]); s=score[tag]; s[1]+=1; s[0]+=int(ok)

for N in (2,3):
    for wv in SEEDS:
        w = wv[:N]
        for e,c in itertools.permutations(range(1,N+1),2):
            qv = w[e-1]/w[c-1]                 # on D: 1 - q w_c/w_e = 0
            if qv in (0,1): continue
            hyp = (lambda M: e in M) if e<c else (lambda M: e in M and c in M)
            name = "e<c" if e<c else "e>c"
            for M in subsets(N):
                if not M: continue
                for I in subsets(N):
                    if len(I)>len(M): continue
                    try:
                        val = sp.simplify(F(N,len(M),len(I),[w[m-1] for m in M],w,I,
                                            sp.Integer(1),sp.Integer(0)).subs({q:qv}))
                    except ZeroDivisionError:
                        continue
                    if val is sp.nan or val.has(sp.zoo): continue
                    if hyp(M): rec(f"N={N} {name}  LEMMA: F vanishes on D", val==0)
                    else:      rec(f"N={N} {name}  PLANT (hyp fails) MUST-NOT-VANISH", val!=0)

print("\n=== F-vanishing lemma, by direct evaluation on D ===")
bad=0
for k in sorted(score):
    g,n=score[k]
    if g!=n: bad+=1
    print(f"  {g:5d}/{n:5d}  {k}" + ("   <<<< DEFECT" if g!=n else ""))
print(f"\n{'ALL GREEN' if bad==0 else str(bad)+' DEFECTIVE LINES'}")
