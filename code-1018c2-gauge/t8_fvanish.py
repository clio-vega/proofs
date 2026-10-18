"""Lemma (F-vanishing at a q-shifted A-type divisor).

D = {w_e = q w_c} = {1 - q w_c/w_e = 0}, e != c.  Claim:
   e < c:  v_D(F_{N,|M|,r}(w_M|w|I)) >= 1  whenever  e in M
   e > c:  v_D(F) >= 1                     whenever  {e,c} subset of M
Proof is a case split on the slot carrying w_e (and, for e>c, also w_c);
see the paper.  PLANTS: the same test with e not in M (resp. {e,c} not subset M)
must FAIL -- if it passed, the lemma would be carrying no information.
"""
import sympy as sp, itertools, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1018c2-gauge')
from q484_engine import q, wvars, F

t = sp.symbols('t')
def subsets(N):
    return [tuple(s) for k in range(N+1) for s in itertools.combinations(range(1,N+1),k)]

def ordt(expr):
    e = sp.cancel(sp.together(sp.expand(expr)))
    if e == 0: return None
    num, den = sp.fraction(e)
    def m(p):
        k = 0
        while True:
            P = sp.Poly(sp.expand(p), t)
            if P.total_degree() == 0: return k
            qq, rr = sp.div(P, sp.Poly(t, t))
            if sp.simplify(rr.as_expr()) == 0: p = qq.as_expr(); k += 1
            else: return k
    return m(num) - m(den)

score = {}
def rec(tag, ok):
    score.setdefault(tag, [0,0]); s = score[tag]; s[1]+=1; s[0]+=int(ok)

for N in (2,3):
    w = wvars(N)
    for e, c in itertools.permutations(range(1, N+1), 2):
        sub = {q: (1-t)*w[e-1]/w[c-1]}       # q = w_e/w_c  on D, t the uniformiser
        hyp = (lambda M: e in M) if e < c else (lambda M: e in M and c in M)
        name = f"e<c" if e < c else f"e>c"
        for M in subsets(N):
            for I in subsets(N):
                if len(I) > len(M): continue
                o = ordt(F(N, len(M), len(I), [w[m-1] for m in M], w, I,
                           sp.Integer(1), sp.Integer(0)).subs(sub))
                if o is None: continue                # identically zero: no information
                if hyp(M):
                    rec(f"N={N} {name}  LEMMA v(F)>=1", o >= 1)
                else:
                    rec(f"N={N} {name}  PLANT (hyp fails) MUST-FAIL", o < 1)

print("\n=== F-vanishing lemma ===")
bad = 0
for k in sorted(score):
    g, n = score[k]
    if g != n: bad += 1
    print(f"  {g:4d}/{n:4d}  {k}" + ("   <<<< DEFECT" if g != n else ""))
print(f"\n{'ALL GREEN' if bad==0 else str(bad)+' DEFECTIVE'}")
