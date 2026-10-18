"""Fast census: exact RATIONAL evaluation at several random seeds.

Same claims as t2 but with q and the free values specialised to exact rationals,
so every comparison is exact arithmetic in Q and finishes.  Three independent
seeds; a claim is scored green only if all seeds agree.  Plants as in t2.
"""
import sympy as sp, itertools, random
from gauge_check import AtilU, FU, EU, pred_Atil_ratio, pred_F_ratio, pred_E_ratio, Q1, Q2, P2
from q484_engine import q, wvars

SEEDS = [(sp.Rational(2,7),  [sp.Rational(3,5), sp.Rational(5,3), sp.Rational(7,2), sp.Rational(11,4)]),
         (sp.Rational(5,11), [sp.Rational(4,9), sp.Rational(9,4), sp.Rational(2,13), sp.Rational(13,6)]),
         (sp.Rational(3,13), [sp.Rational(6,5), sp.Rational(5,6), sp.Rational(8,3), sp.Rational(17,5)])]

score = {}
def rec(tag, ok):
    score.setdefault(tag,[0,0]); s=score[tag]; s[1]+=1; s[0]+=int(ok)

def Dsets(xv, ev):
    return [("d=1 fixed-pt", [ev]),
            ("d=2 free orbit {x,1/x}", [xv, 1/xv]),
            ("d=2 two fixed pts {e,e}", [ev, ev]),
            ("d=3 {x,1/x,-1}", [xv, 1/xv, sp.Integer(-1)]),
            ("d=4 {x,1/x,x,1/x}", [xv, 1/xv, xv, 1/xv]),
            ("d=4 {x,1/x,y,1/y}", [xv, 1/xv, sp.Rational(7,3), sp.Rational(3,7)])]

for si,(qv, vals) in enumerate(SEEDS):
    xv, r1, r2, w2v = vals
    for N in (1,2):
        wl = [vals[3], sp.Rational(5,7)][:N]
        for ev in (sp.Integer(1), sp.Integer(-1)):
            for dn, D in Dsets(xv, ev):
                for Rn, R in [("R=[]", []), ("R=[r1]", [r1]), ("R=[r1,r2]", [r1, r2])]:
                    for k in range(N+1):
                        for I in itertools.combinations(range(1,N+1), k):
                            if len(R) < len(I): continue
                            if len(R)+len(D) > 5: continue
                            sub = {q: qv}
                            aU = sp.nsimplify(AtilU(N, [v for v in R]+list(D), wl, I, 'C').subs(sub))
                            aR = sp.nsimplify(AtilU(N, R, wl, I, 'C').subs(sub))
                            if sp.simplify(aR) == 0:
                                rec(f"{dn} [vacuous stratum]", sp.simplify(aU)==0); continue
                            pred = pred_Atil_ratio(D, 2).subs(sub)
                            rec(f"A {dn}  Atil-ratio", sp.simplify(aU/aR - pred)==0)
                            # F and E ratios SEPARATELY (a quotient can be 1 with both wrong)
                            fU, fR = FU(N,[v for v in R]+list(D),wl,I,'C').subs(sub), FU(N,R,wl,I,'C').subs(sub)
                            if sp.simplify(fR)!=0:
                                rec(f"A {dn}  F-ratio", sp.simplify(fU/fR - pred_F_ratio(N,R,D,wl).subs(sub))==0)
                            eU, eR = EU(N,[v for v in R]+list(D),wl,I,'C').subs(sub), EU(N,R,wl,I,'C').subs(sub)
                            rec(f"A {dn}  E-ratio", sp.simplify(eU/eR - pred_E_ratio(N,R,D,wl,'C').subs(sub))==0)
    # --------- plants: D NOT iota-closed, or prediction at the wrong h
    N=2; wl=[vals[3], sp.Rational(5,7)]
    PL = [("P1 D=[x,y], xy!=1", [xv, sp.Rational(7,3)], 2),
          ("P2 D=[2], 2^2!=1",  [sp.Integer(2)], 2),
          ("P3 D=[x,1/x] at h=1 (wrong h)", [xv, 1/xv], 1),
          ("P4 D=[x,x], not iota-closed", [xv, xv], 2)]
    for pn, D, hh in PL:
        for R in ([], [r1]):
            for I in ((), (1,)):
                if len(R) < len(I): continue
                aU = AtilU(N, list(R)+list(D), wl, I, 'C').subs({q:qv})
                aR = AtilU(N, R, wl, I, 'C').subs({q:qv})
                if sp.simplify(aR)==0: continue
                pred = pred_Atil_ratio(D, hh).subs({q:qv})
                rec(f"C {pn} MUST-FAIL", sp.simplify(aU/aR - pred)!=0)

# --------- identities, symbolic (cheap)
x,y = sp.symbols('x y')
def idents(D, must):
    d=len(D); py=sp.Integer(1)
    for v in D: py*=v
    den=sp.Integer(1)
    for v in D: den*=(1-q*v**2)
    rec(("D P2=(prod y)^(d-1) Q2" if must else "D P2 ident MUST-FAIL"),
        (sp.simplify(P2(D)-py**(d-1)*Q2(D))==0)==must)
    rec(("D Q2^2/Q1=1/prod(1-q y^2)" if must else "D Q2^2/Q1 MUST-FAIL"),
        (sp.simplify(Q2(D)**2/Q1(D)-1/den)==0)==must)
for D in ([sp.Integer(1)],[sp.Integer(-1)],[x,1/x],[sp.Integer(1),sp.Integer(1)],
          [sp.Integer(-1),sp.Integer(-1)],[x,1/x,sp.Integer(-1)],[x,1/x,x,1/x],[x,1/x,y,1/y]):
    idents(D, True)
for D in ([x,y],[x,x],[sp.Integer(2)]):
    idents(D, False)

print("\n================ FAST CENSUS ================")
bad=0
for k in sorted(score):
    g,n = score[k]
    if g!=n: bad+=1
    print(f"  {g:5d}/{n:5d}  {k}" + ("   <<<< DEFECT" if g!=n else ""))
print(f"\n{'ALL GREEN' if bad==0 else str(bad)+' DEFECTIVE LINES'}")
