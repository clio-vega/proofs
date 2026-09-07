import sys, sympy as sp
sys.path.insert(0,'/home/clio/projects/scratch/q92')
import engineA as A, engineB as B
t = sp.Symbol('t')

def split(lam, e, f):
    """separate T1 (one-bead) and T2 (two-bead) contributions of engine B"""
    M = B.maya(lam); out1={}; out2={}
    for b in range(-B.W+1, B.W):
        if b+e+f >= B.W: continue
        d = B.occ(M,b+f)-B.occ(M,b+e)
        if d==0: continue
        st = B.ann((sp.Integer(1),M), b)
        if st is None: continue
        st = B.cre(st, b+e+f)
        if st is None: continue
        s,Mp = st
        c = -(1+t)/t * s * (-t)**B.Ncount(M,b,b+e+f) * d
        k=B.unmaya(Mp); out1[k]=sp.expand(out1.get(k,0)+c)
    for b in range(-B.W+1,B.W):
        for c_ in range(-B.W+1,B.W):
            if b==c_ or b+e>=B.W or c_+f>=B.W: continue
            k = (1 if c_<b+e<c_+f else 0)-(1 if c_<b<c_+f else 0)
            if k==0: continue
            st=B.ann((sp.Integer(1),M),c_)
            if st is None: continue
            st=B.ann(st,b)
            if st is None: continue
            st=B.cre(st,c_+f)
            if st is None: continue
            st=B.cre(st,b+e)
            if st is None: continue
            s,Mp=st
            co = -(t**2-1)/t*k*s*(-t)**B.Ncount(M,b,b+e)*(-t)**B.Ncount(M,c_,c_+f)
            kk=B.unmaya(Mp); out2[kk]=sp.expand(out2.get(kk,0)+co)
    f1={k:v for k,v in out1.items() if sp.expand(v)!=0}
    f2={k:v for k,v in out2.items() if sp.expand(v)!=0}
    return f1,f2

# beads moved, from engine A side: |M(mu) sym-diff M(lam)|/2
def nmoved(lam,mu):
    return len(B.maya(lam) ^ B.maya(mu))//2

print("=== Psi_{1,f} should vanish identically (predicted: k=0 whenever Q!=0) ===")
for f in range(2,7):
    tot=set()
    for n in range(0,7):
        for lam in A.partitions(n):
            _,p2 = split(lam,1,f)
            tot |= set(p2)
    print(f"  Psi_(1,{f}) support: {sorted(tot) if tot else 'EMPTY'}")

print("\n=== two-bead sector of [R_e,R_f] read off ENGINE A (independent) ===")
for (e,f) in [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(2,5)]:
    twob=[]
    for n in range(0,7):
        for lam in A.partitions(n):
            for mu,c in A.commutator(lam,e,f).items():
                if nmoved(lam,mu)==2: twob.append((lam,mu,c))
    print(f"  (e,f)=({e},{f}): {len(twob)} nonzero 2-bead matrix elements", 
          ("  e.g. "+str(twob[0])) if twob else "")
