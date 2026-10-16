"""The residual hypothesis of the converse.

LEMMA (nonvanH).  Let a<b, H = H_ab = {w_a w_b = 1}.  For M >= I with
{a,b} cap M = empty,   ord_H Atilde_{M,I} = 0,  i.e. Atilde_{M,I}|_H != 0.

This is exactly the H-LOCALISATION of the converse of lem:tri (M>=I =>
A_{M,I} != 0), which Q484 measured 84/84.  Report as a set, with the two
domains named:
  D0 : {a,b} cap M = empty          -- what the converse of the theorem needs
  D1 : {a,b} not subset M (so |{a,b} cap M| <= 1)  -- the wider claim

CONTROL (positive): also report, over the SAME loop, the pairs with M >= I and
{a,b} subset M.  There lem:A predicts ord_H >= 1 whenever M\{a,b} not>= I, and
the count of ord_H = 0 cases must be exactly the M\{a,b} >= I ones.  A scan
that cannot distinguish the two domains is decoration.
"""
import sys; sys.path.insert(0,'.')
import sympy as sp, time
from collections import Counter
import q484_engine as E
q = E.q

def ordH(expr, w, a, b, maxk=4):
    """ord of expr along H_ab (expr a Laurent polynomial => ord >= 0)."""
    expr = sp.cancel(sp.together(expr))
    if expr == 0: return sp.oo
    s = sp.Symbol('s')
    num,den = sp.fraction(expr)
    num_s = sp.expand(num.subs(w[b-1], s/w[a-1]))
    den_s = sp.expand(den.subs(w[b-1], s/w[a-1]))
    pn = sp.Poly(num_s, s); pd = sp.Poly(den_s, s)
    def mult(p):
        k=0; pp=p
        while k<=maxk:
            if pp.eval(1) != 0: return k
            pp = sp.Poly(sp.cancel(pp.as_expr()/(s-1)), s); k+=1
        return k
    return mult(pn)-mult(pd)

def run(N, X='C'):
    w = E.wvars(N); S = E.subsets(N)
    t0=time.time()
    A = {(M,I): E.Atilde(N,w,M,I,X) for M in S for I in S}
    print(f'=== nonvanH  N={N} {X}  (build {time.time()-t0:.0f}s) ===',flush=True)
    stats = {'D0':Counter(), 'D1only':Counter(), 'IN':Counter()}
    viol = []
    ctrl_ok=0; ctrl_bad=0; ctrl_viol=[]
    for a in range(1,N+1):
      for b in range(a+1,N+1):
        for M in S:
          for I in S:
            if not E.geq(M,I): continue
            inter = len({a,b} & set(M))
            o = ordH(A[(M,I)], w, a, b)
            if inter == 2:
                Mp = tuple(x for x in M if x not in (a,b))
                pred0 = E.geq(Mp,I)        # lem:A: ord>=1 unless M' >= I
                stats['IN'][o]+=1
                if (o == 0) == pred0: ctrl_ok+=1
                else:
                    ctrl_bad+=1
                    if len(ctrl_viol)<8: ctrl_viol.append((a,b,M,I,o,pred0))
            elif inter == 0:
                stats['D0'][o]+=1
                if o != 0 and len(viol)<10: viol.append(('D0',a,b,M,I,o))
            else:
                stats['D1only'][o]+=1
                if o != 0 and len(viol)<10: viol.append(('D1only',a,b,M,I,o))
    for k in ('D0','D1only','IN'):
        print(f'  {k:7s} ord_H distribution: {dict(stats[k])}  (n={sum(stats[k].values())})')
    n0 = sum(stats['D0'].values()); z0 = stats['D0'][0]
    n1 = sum(stats['D1only'].values()); z1 = stats['D1only'][0]
    print(f'\n  LEMMA nonvanH on D0 ({{a,b}} cap M = empty): {z0}/{n0} have ord_H = 0')
    print(f'  LEMMA nonvanH on D1only (|cap|=1)           : {z1}/{n1} have ord_H = 0')
    print(f'  violations: {viol if viol else "none"}')
    print(f'\n  POSITIVE CONTROL (|cap|=2, lem:A sharpness): {ctrl_ok} agree, {ctrl_bad} disagree')
    print(f'    -- the control DISCRIMINATES: ord_H distribution above has both 0 and >=1')
    for z in ctrl_viol: print('      ctrl mismatch a,b,M,I,ord,pred(M\'>=I) =',z)

if __name__=='__main__':
    run(int(sys.argv[1]), sys.argv[2] if len(sys.argv)>2 else 'C')
