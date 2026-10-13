"""PRE-REGISTERED SCREEN (C) + the main verification.  Domains printed with every verdict."""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1013c2-q442')
from poly import *
from q442 import *

def nz(p): return len(p)
def eq(p, q): return p == q
def nonneg(p): return all(c > 0 for c in p.values())

DOMAINS = [(1,0,2),(1,0,3),(1,0,4),(1,1,2),(1,1,3),(1,1,4),(1,1,5),(1,2,4),(1,2,5),
           (1,3,4),(2,0,2),(2,0,3),(2,0,4),(2,1,3),(2,1,4),(2,2,3),(2,2,4),(3,0,3),(3,1,3)]

print("="*100)
print("Q442 SCREEN.  M=k+ell+1, N=2M, w=2ell+2.  Engines: A=det, B=Warnaar Laurent, C=path DP")
print("="*100)
hdr = f"{'(k,ell,n)':>12} {'w':>3} {'N':>3} | {'A=B':>5} {'A=A^T':>6} {'A=C:Ainf':>9} | " \
      f"{'[y0]=[k0]?':>10} {'#mon diff':>9} {'diff>=0':>8} {'diff=ow+Binf':>13} | " \
      f"{'t=1 gate':>9} {'#tdeg':>6}"
print(hdr); print("-"*len(hdr))

rows = []
for (k, ell, n) in DOMAINS:
    w = 2*ell + 2; M = k+ell+1; N = 2*M
    A  = engineA(k, n)
    AT = engineA_transposed(k, n)
    B, G1, tdegs = engineB(k, ell, n)
    C  = engineC_sums(k, n, w)
    y0 = A
    k0 = C['kappa0']
    diff = sub(y0, k0)
    pred = add(C['overwide'], C['Binf'])
    rows.append(dict(k=k,ell=ell,n=n,w=w,
        AB=eq(A,B), AAT=eq(A,AT), ACinf=eq(A,C['Ainf']),
        sector=eq(y0,k0), ndiff=nz(diff), diffpos=(nonneg(diff) if diff else None),
        diffpred=eq(diff,pred),
        t1=eq(G1, sub(C['A'], C['B'])), ntdeg=len(tdegs),
        AeqA=eq(C['Ainf'], C['A'])))
    r = rows[-1]
    print(f"{str((k,ell,n)):>12} {w:>3} {N:>3} | {str(r['AB']):>5} {str(r['AAT']):>6} "
          f"{str(r['ACinf']):>9} | {str(r['sector']):>10} {r['ndiff']:>9} "
          f"{str(r['diffpos']):>8} {str(r['diffpred']):>13} | {str(r['t1']):>9} {r['ntdeg']:>6}")

print("-"*len(hdr))
tot = len(rows)
print(f"DOMAIN OF EVERY VERDICT BELOW: the {tot} triples (k,ell,n) listed above, "
      f"k in {{1,2,3}}, ell in {{0,1,2,3}}, n <= 5.")
print()
print(f"  A == B  (det engine == Warnaar [t^0] engine)          : {sum(r['AB'] for r in rows)}/{tot}")
print(f"  A == A^T                                              : {sum(r['AAT'] for r in rows)}/{tot}")
print(f"  A == C.Ainf  (= sum_alpha A^inf_alpha x^alpha)         : {sum(r['ACinf'] for r in rows)}/{tot}")
print(f"  t=1 gate: G(1)*(x..)^k == C.A - C.B                   : {sum(r['t1'] for r in rows)}/{tot}")
print()
print(f"  SCREEN (C): [y=0] == [kappa=0]                         : "
      f"{sum(r['sector'] for r in rows)}/{tot}   <-- FIRES where False")
print(f"  of the {tot}, number where the two DIFFER             : {sum(not r['sector'] for r in rows)}")
print(f"  difference has all coefficients > 0 wherever nonzero  : "
      f"{sum(1 for r in rows if r['diffpos'] is True)}/{sum(not r['sector'] for r in rows)}")
print(f"  difference == overwide-even class + B^inf class        : {sum(r['diffpred'] for r in rows)}/{tot}")
print()
print(f"  NON-VACUITY of the t-grading (#distinct t-degrees in G): "
      f"min {min(r['ntdeg'] for r in rows)}, max {max(r['ntdeg'] for r in rows)}; "
      f"rows with only one t-degree: {sum(1 for r in rows if r['ntdeg']<2)}")
print(f"  NON-VACUITY of the screen: rows where [y=0]==[kappa=0] (screen blind): "
      f"{[ (r['k'],r['ell'],r['n']) for r in rows if r['sector']]}")
print()
print("  ARM 1 control  A^inf == A (width bound inactive):")
for r in rows:
    flag = "n<=w" if r['n'] <= r['w'] else "n> w"
    print(f"      (k,ell,n)=({r['k']},{r['ell']},{r['n']})  w={r['w']}  {flag}  "
          f"A^inf==A: {r['AeqA']}")
