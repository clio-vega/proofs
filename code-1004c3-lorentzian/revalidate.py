"""RE-VALIDATION after fixing count_pos_eigen_exact.  The defective version was the
independent validator of l3-det-reduction, so that validation must be redone -- and it
must now be made to see repeated eigenvalues, the case it was blind to."""
import random, itertools
from lor import *
random.seed(20261004)

print("=== V1'. l3-det-reduction vs the CORRECTED exact count ===")
n=dis=0
for _ in range(600):
    r=lambda: random.randint(0,12)
    a,b,c0,d,e,f=r(),r(),r(),r(),r(),r()
    M=[[a,b,c0],[b,d,e],[c0,e,f]]
    n+=1
    if at_most_one_pos_3x3_detcrit(M)!=(count_pos_eigen_exact(M)<=1): dis+=1; print("  DISAGREE",M)
print(f"  random nonnegative symmetric 3x3 enumerated {n}, disagreements {dis}")

print("=== V1''. EXHAUSTIVE over small entries, and separately over matrices with")
print("===       REPEATED eigenvalues -- the class the defective validator was blind to ===")
n=dis=0; rep=repdis=0
for a,b,c0,d,e,f in itertools.product(range(0,5),repeat=6):
    M=[[a,b,c0],[b,d,e],[c0,e,f]]
    n+=1
    ex=count_pos_eigen_exact(M)
    if at_most_one_pos_3x3_detcrit(M)!=(ex<=1): dis+=1; print("  DISAGREE",M)
    import sympy as sp
    if len(sp.Matrix(M).charpoly(sp.Symbol('l')).as_poly().real_roots())>len(set(sp.Matrix(M).charpoly(sp.Symbol('l')).as_poly().real_roots())):
        rep+=1
        if at_most_one_pos_3x3_detcrit(M)!=(ex<=1): repdis+=1
print(f"  exhaustive entries 0..4 enumerated {n}, disagreements {dis}")
print(f"  of these, matrices with a REPEATED eigenvalue: {rep}, disagreements among them {repdis}")

print("=== Lemma minor, re-run with the corrected counter ===")
tot=ok=viol=0; wit=[]
for n_ in (2,3,4,5):
    for _ in range(500):
        A=[[random.randint(-6,6) for _ in range(n_)] for _ in range(n_)]
        M=[[A[i][j]+A[j][i] for j in range(n_)] for i in range(n_)]
        if any(M[i][i]<0 for i in range(n_)): continue
        if count_pos_eigen_exact(M)>1: continue
        tot+=1
        bad=[(i,j) for i in range(n_) for j in range(n_) if M[i][i]*M[j][j]-M[i][j]**2>0]
        if bad: viol+=1; wit.append(M)
        else: ok+=1
# plus a deliberate sweep over SMALL integer matrices including repeated-eigenvalue ones
tot2=viol2=0
for a,b,c0,d,e,f in itertools.product(range(0,4),repeat=6):
    for sgn in (1,-1):
        M=[[a,sgn*b,c0],[sgn*b,d,e],[c0,e,f]]
        if any(M[i][i]<0 for i in range(3)): continue
        if count_pos_eigen_exact(M)>1: continue
        tot2+=1
        if any(M[i][i]*M[j][j]-M[i][j]**2>0 for i in range(3) for j in range(3)): viol2+=1
print(f"  random, n=2..5: hypotheses satisfied by {tot}, conclusion held {ok}, violations {viol}")
print(f"  exhaustive 3x3, entries 0..3 (both signs off-diagonal): satisfied by {tot2}, violations {viol2}")
if wit: print("  WITNESSES:", wit[:3])
