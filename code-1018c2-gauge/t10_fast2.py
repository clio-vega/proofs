"""Fast census, with a GENERICITY GUARD.

t9 reported 11 defective lines.  Every one was a degenerate seed, not a defect:
  - seed 1 had r1 = 5/3 and D = [3/5, 5/3], so a value of R EQUALLED a value of
    D; two u-arguments coincide and the (v_j - v_k) denominators of lines 2a/3a
    blow up (nan), which is outside the theorem's "generic point" hypothesis;
  - seed 1 also had q*w_1 = (2/7)(7/2) = 1, so prod_k(1 - q y w_k) = 0 and the
    predicted F-ratio was 0/0.
The symbolic smoke test was green on the same claims because symbols are generic
by construction.  So: validate the seed, and REPORT the skipped instances rather
than dropping them silently.
"""
import sympy as sp, itertools
from gauge_check import AtilU, FU, EU, pred_Atil_ratio, pred_F_ratio, pred_E_ratio, Q1, Q2, P2
from q484_engine import q, wvars

def generic(qv, vals, wl):
    """no value collisions, no inverse collisions, no q*monomial = 1."""
    pool = list(vals) + list(wl)
    for i,a in enumerate(pool):
        for b in pool[i+1:]:
            if a == b or a*b == 1: return False
        if a**2 == 1 and a not in (1,-1): return False
    for a in pool:
        if qv*a == 1: return False
        for b in pool:
            if qv*a*b == 1 or qv*a/b == 1: return False
    return True

SEEDS = [(sp.Rational(2,7),  [sp.Rational(3,5), sp.Rational(7,11), sp.Rational(13,17)], [sp.Rational(19,23), sp.Rational(29,31)]),
         (sp.Rational(5,13), [sp.Rational(4,9), sp.Rational(11,7),  sp.Rational(23,19)], [sp.Rational(37,41), sp.Rational(43,47)]),
         (sp.Rational(3,11), [sp.Rational(6,5), sp.Rational(17,13), sp.Rational(31,29)], [sp.Rational(53,59), sp.Rational(61,67)])]

score = {}; skipped = 0
def rec(tag, ok):
    score.setdefault(tag,[0,0]); s=score[tag]; s[1]+=1; s[0]+=int(ok)

for qv, vals, wfull in SEEDS:
    xv, r1, r2 = vals
    assert generic(qv, vals, wfull), "seed itself not generic"
    for N in (1,2):
        wl = wfull[:N]
        for ev in (sp.Integer(1), sp.Integer(-1)):
            DS = [("d=1 fixed-pt", [ev]),
                  ("d=2 free orbit {x,1/x}", [xv, 1/xv]),
                  ("d=2 two fixed pts {e,e}", [ev, ev]),
                  ("d=3 {x,1/x,eps}", [xv, 1/xv, ev]),
                  ("d=4 {x,1/x,x,1/x}", [xv, 1/xv, xv, 1/xv])]
            for dn, D in DS:
                for R in ([], [r1], [r1, r2]):
                    if len(R)+len(D) > 5: continue
                    # genericity of THIS instance: R-values and D-values must not collide
                    bad = False
                    for a in R:
                        for b in list(D)+list(wl):
                            if a==b or a*b==1 or qv*a*b==1 or qv*a/b==1 or qv*b/a==1: bad=True
                    for b in D:
                        for c in wl:
                            if qv*b*c==1: bad=True
                    if bad:
                        skipped += 1; continue
                    for k in range(N+1):
                        for I in itertools.combinations(range(1,N+1),k):
                            if len(R) < len(I): continue
                            U = list(R)+list(D)
                            aR = AtilU(N, R, wl, I,'C').subs({q:qv})
                            aU = AtilU(N, U, wl, I,'C').subs({q:qv})
                            if sp.simplify(aR)==0:
                                rec(f"{dn} [vacuous stratum]", sp.simplify(aU)==0); continue
                            pred = sp.simplify(pred_Atil_ratio(D,2).subs({q:qv}))
                            rec(f"A {dn}  Atil-ratio", sp.simplify(sp.simplify(aU/aR)-pred)==0)
                            fU,fR = FU(N,U,wl,I,'C').subs({q:qv}), FU(N,R,wl,I,'C').subs({q:qv})
                            if sp.simplify(fR)!=0:
                                rec(f"A {dn}  F-ratio",
                                    sp.simplify(sp.simplify(fU/fR)-sp.simplify(pred_F_ratio(N,R,D,wl).subs({q:qv})))==0)
                            eU,eR = EU(N,U,wl,I,'C').subs({q:qv}), EU(N,R,wl,I,'C').subs({q:qv})
                            rec(f"A {dn}  E-ratio",
                                sp.simplify(sp.simplify(eU/eR)-sp.simplify(pred_E_ratio(N,R,D,wl,'C').subs({q:qv})))==0)
    # plants, same generic seed
    N=2; wl=wfull
    PL = [("P1 D=[x,y], xy!=1", [xv, r2], 2),
          ("P2 D=[2], 2^2!=1",  [sp.Integer(2)], 2),
          ("P3 D=[x,1/x] scored at h=1", [xv, 1/xv], 1),
          ("P4 D=[x,x], not iota-closed", [xv, xv], 2)]
    for pn, D, hh in PL:
        for R in ([], [r1]):
            for I in ((), (1,)):
                if len(R)<len(I): continue
                aR = AtilU(N,R,wl,I,'C').subs({q:qv})
                if sp.simplify(aR)==0: continue
                aU = AtilU(N,list(R)+list(D),wl,I,'C').subs({q:qv})
                pred = sp.simplify(pred_Atil_ratio(D,hh).subs({q:qv}))
                rec(f"C {pn} MUST-FAIL", sp.simplify(sp.simplify(aU/aR)-pred)!=0)

x,y = sp.symbols('x y')
def idents(D, must, admissible=True):
    d=len(D); py=sp.Integer(1)
    for v in D: py*=v
    den=sp.Integer(1)
    for v in D: den*=(1-q*v**2)
    tag1 = "D P2=(prod y)^(d-1) Q2" if must else "D P2 ident MUST-FAIL"
    tag2 = "D Q2^2/Q1=1/prod(1-q y^2)" if must else "D Q2^2/Q1 MUST-FAIL"
    if admissible:
        rec(tag1, (sp.simplify(P2(D)-py**(d-1)*Q2(D))==0)==must)
    rec(tag2, (sp.simplify(Q2(D)**2/Q1(D)-1/den)==0)==must)
for D in ([sp.Integer(1)],[sp.Integer(-1)],[x,1/x],[sp.Integer(1),sp.Integer(1)],
          [sp.Integer(-1),sp.Integer(-1)],[x,1/x,sp.Integer(-1)],[x,1/x,x,1/x],[x,1/x,y,1/y]):
    idents(D, True)
# plants.  NOTE: D=[2] is INADMISSIBLE for the P2 identity -- at d=1 both sides are
# empty products and the identity reads 1=1, so it cannot fail.  t9 scored it and got
# 2/3; that was the plant being inadmissible, not the identity being wrong.
idents([x,y],  False)
idents([x,x],  False)
idents([sp.Integer(2)], False, admissible=False)

print("\n============ FAST CENSUS (generic seeds) ============")
bad=0
for k in sorted(score):
    g,n = score[k]
    if g!=n: bad+=1
    print(f"  {g:5d}/{n:5d}  {k}" + ("   <<<< DEFECT" if g!=n else ""))
print(f"\nskipped as non-generic (reported, not dropped): {skipped} (R,D,w) combinations")
print(f"{'ALL GREEN' if bad==0 else str(bad)+' DEFECTIVE LINES'}")
