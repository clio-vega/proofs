"""Instrument validation: every instrument must reproduce a KNOWN PRESENT value
before any control or sweep is believed."""
import random, itertools, sys
from lor import *
from lor import _char_poly_sym3

random.seed(20261004)
print("=== V1. det-criterion vs exact Sturm count, nonnegative symmetric 3x3 ===")
n_tested = n_dis = 0
for _ in range(600):
    r = lambda: random.randint(0, 12)
    a,b,c0,d,e,f = r(),r(),r(),r(),r(),r()
    M = [[a,b,c0],[b,d,e],[c0,e,f]]
    crit = at_most_one_pos_3x3_detcrit(M)
    exact = count_pos_eigen_exact(M) <= 1
    n_tested += 1
    if crit != exact:
        n_dis += 1; print("  DISAGREE", M, crit, exact)
print(f"  matrices enumerated: {n_tested}, disagreements: {n_dis}")

print("=== V2. the rlc-implies-l3 witness: 2x2 minors all <=0 but det<0 ===")
# a symmetric 3x3, nonneg entries, all three 2x2 principal minors <=0, det = -16
found = None
for a,b,c0,d,e,f in itertools.product(range(0,7),repeat=6):
    M=[[a,b,c0],[b,d,e],[c0,e,f]]
    m01=a*d-b*b; m02=a*f-c0*c0; m12=d*f-e*e
    _,_,det=_char_poly_sym3(M)
    if m01<=0 and m02<=0 and m12<=0 and det==-16:
        found=M; break
print("  witness:", found)
if found:
    print("  all 2x2 principal minors <=0 :", True)
    print("  det =", _char_poly_sym3(found)[2])
    print("  det-criterion says at-most-one-pos:", at_most_one_pos_3x3_detcrit(found), "(must be False)")
    print("  exact positive eigenvalue count  :", count_pos_eigen_exact(found), "(must be >=2)")

print("=== V3. is_Mconvex on known present / known absent sets ===")
S_box = [a for a in itertools.product(range(0,4),repeat=3) if sum(a)==4 and a[2]<=2]
print("  box slice {|a|=4, a_2<=2}: size", len(S_box), "M-convex:", is_Mconvex(S_box), "(must be True)")
S_bad = [(2,0,0),(0,2,0),(0,0,2)]
print("  {2e_i} (no exchange):      size", len(S_bad), "M-convex:", is_Mconvex(S_bad), "(must be False)")
S_bad2 = [(2,0,0),(0,2,0),(1,1,0),(0,0,2)]
print("  that + (1,1,0):            size", len(S_bad2),"M-convex:", is_Mconvex(S_bad2),"(must be False)")

print("=== V4. k_sequence vs direct expansion of F(z)=sum_w prod [w_i]_z ===")
def F_direct(L,H,T):
    m=len(L)
    out={}
    for w in itertools.product(*[range(L[i],H[i]+1) for i in range(m)]):
        if sum(w)!=T: continue
        p={0:1}
        for wi in w:
            q={}
            for k,v in p.items():
                for t in range(wi): q[k+t]=q.get(k+t,0)+v
            p=q
        for k,v in p.items(): out[k]=out.get(k,0)+v
    if not out: return []
    return [out.get(a,0) for a in range(max(out)+1)]
agree=tot=0
for m in (1,2,3):
    for L in itertools.product(range(1,4),repeat=m):
        for g in itertools.product(range(0,3),repeat=m):
            H=tuple(L[i]+g[i] for i in range(m))
            for T in range(sum(L), sum(H)+1):
                A=F_direct(L,H,T)                  # in w-language, w_i>=1
                B=k_sequence([x-1 for x in L],[x-1 for x in H], T-m)  # in u-language
                A=A+[0]*(max(0,len(B)-len(A))); B=B+[0]*(max(0,len(A)-len(B)))
                tot+=1
                if A==B: agree+=1
                elif tot<2000 and agree<3: print("   MISMATCH",L,H,T,A,B)
print(f"  instances enumerated: {tot}, agreements: {agree}")

print("=== V5. prop:m3 known value  L=U=(2,2,2) in u-language, D=6 ===")
ks=k_sequence([2,2,2],[2,2,2],6)
print("  k =", ks, " (must be [1,3,6,7,6,3,1])")
print("  PF2:", is_logconcave_pf2(ks), "  concave:", all(2*ks[a]>=ks[a-1]+ks[a+1] for a in range(1,len(ks)-1)))
