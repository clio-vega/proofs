"""Q255 ell=3: test conditions (A) RLC and (B) M_ik M_kj >= M_ij M_kk on the
Hessians M_ij = K^c_{beta+e_i+e_j}.  All exact integer arithmetic."""
import sys, itertools
sys.path.insert(0,'/home/clio/projects/proofs/code-q254')
sys.path.insert(0,'/home/clio/projects/proofs/code-2026-09-19')
sys.path.insert(0,'/home/clio/projects/proofs/code-q255')
import cyl as C
from winding import weights_counted, winds
from lorentzian import is_M_convex, compositions
from falsify import cyl_shapes, n_positive_eigs_exact

def hessians(K, ell, d):
    for beta in compositions(d-2, ell):
        M=[[K.get(tuple(g),0) for g in
            [tuple(x+ (1 if t==i else 0)+(1 if t==j else 0) for t,x in enumerate(beta))
             for j in range(ell)]] for i in range(ell)]
        yield beta, M

def condA(M, ell):           # M_ij^2 >= M_ii M_jj      (i != j)
    return [(i,j) for i in range(ell) for j in range(ell) if i!=j
            and M[i][j]**2 < M[i][i]*M[j][j]]

def condB(M, ell):           # M_ik M_kj >= M_ij M_kk    (i,j,k distinct)
    return [(i,j,k) for i,j,k in itertools.permutations(range(ell),3)
            if M[i][k]*M[k][j] < M[i][j]*M[k][k]]

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
           -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
           +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

if __name__=="__main__":
    ell=3
    tot=hess=wind_h=0
    failA=failB=failDet=failL3=0
    exA=[];exB=[];exD=[]
    zero_det=0
    for n in range(2,7):
      for m in range(1,n+1):
        for mu in cyl_shapes(n,m):
          for lam in itertools.product(*[range(mu[i],mu[i]+n+1) for i in range(m)]):
            if not C.is_shape(lam,n,m) or not C.contains(lam,mu): continue
            d=C.size(mu,lam)
            if d<2 or d>9: continue
            K=weights_counted(mu,lam,n,m,ell)
            K={a:c for a,c in K.items() if c>0}
            if not K: continue
            w,_,_=winds(mu,lam,n,m,ell)
            tot+=1
            for beta,M in hessians(K,ell,d):
                if all(M[i][j]==0 for i in range(3) for j in range(3)): continue
                hess+=1; wind_h+= 1 if w else 0
                a=condA(M,3); b=condB(M,3); dt=det3(M)
                npos=n_positive_eigs_exact(M)
                if a: failA+=1;  exA.append((n,m,mu,lam,d,beta,M,a[0]))
                if b: failB+=1;  exB.append((n,m,mu,lam,d,beta,M,b[0],npos,dt))
                if dt<0: failDet+=1; exD.append((n,m,mu,lam,d,beta,M,dt))
                if dt==0: zero_det+=1
                if npos>1: failL3+=1
    print("ell=3 instances=%d  nonzero Hessians=%d (winding-shape: %d)"%(tot,hess,wind_h))
    print("  (A) RLC      failures: %d"%failA)
    print("  (B) M_ik M_kj >= M_ij M_kk failures: %d"%failB)
    print("  det M < 0    failures: %d"%failDet)
    print("  (L3) >1 positive eig  : %d"%failL3)
    print("  det M == 0            : %d  (%.1f%%)"%(zero_det,100*zero_det/max(1,hess)))
    for tag,ex in (("A",exA),("B",exB),("D",exD)):
        for e in ex[:4]: print("   %s-fail"%tag, e)
