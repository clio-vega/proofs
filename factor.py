import sympy as sp, sys
sys.path.insert(0, '/home/clio/projects/reviews/2026-09-08-selfreview-code')
from ribbon import border_strips, parts
t, s = sp.symbols('t s')

def conj(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))
def basis(N): return sorted(parts(N))
def Rmat(e,N,u):
    src,tgt=basis(N),basis(N+e); M=sp.zeros(len(tgt),len(src))
    for j,l in enumerate(src):
        for (mu,h) in border_strips(l,e): M[tgt.index(mu),j]+=u**h
    return M
def maya(lam,L):
    lam=list(lam)+[0]*(L-len(lam)); return frozenset(lam[j]-(j+1)+L for j in range(L))

# --- (A) algebraic identity: the 4-term e=f sum factors ---
P,Q,k=sp.symbols('P Q k')
for Pv in range(0,6):
  for Qv in range(0,6):
    for kv in (-2,-1,0,1,2):
      four = t**(Pv-kv)*s**Qv - t**Pv*s**(Qv+kv) + t**(Qv+kv)*s**Pv - t**Qv*s**(Pv-kv)
      fac  = (1-(t*s)**kv)*(t**(Pv-kv)*s**Qv - t**Qv*s**(Pv-kv))
      assert sp.simplify(sp.expand(four-fac))==0, (Pv,Qv,kv)
print("(A) 4-term sum == (1-(ts)^k)(t^{P-k}s^Q - t^Q s^{P-k}) : verified for P,Q in 0..5, k in -2..2")

# --- (B) two-bead sector vanishes on ts=1, for all e,f (not just e=f) ---
print("\n(B) two-bead part of [R_e(t),R_f(s)] restricted to ts=1:")
for (e,f) in [(1,2),(2,2),(2,3),(3,3),(3,4),(2,5),(4,4)]:
    worst=0; nz2=0
    for N in range(0,5):
        C = Rmat(e,N+f,t)*Rmat(f,N,s) - Rmat(f,N+e,s)*Rmat(e,N,t)
        src,tgt = basis(N), basis(N+e+f)
        L=N+e+f+4
        for j,lam in enumerate(src):
            Ml=maya(lam,L)
            for i,mu in enumerate(tgt):
                if C[i,j]==0: continue
                if len(Ml ^ maya(mu,L))==4:            # two-bead
                    nz2+=1
                    val=sp.simplify(sp.expand(C[i,j].subs(s,1/t)))
                    if val!=0: worst+=1
    print(f"   (e,f)=({e},{f}): {nz2} nonzero two-bead entries, "
          f"{worst} survive on ts=1   -> {'VANISHES' if worst==0 else 'DOES NOT VANISH'}")

# --- (C) HONEST LIMIT: is the one-bead part also zero on ts=1? ---
print("\n(C) one-bead part on ts=1 (the control that must MOVE):")
for (e,f) in [(2,3),(3,3),(3,4)]:
    nz1=0; surv=0
    for N in range(0,5):
        C = Rmat(e,N+f,t)*Rmat(f,N,s) - Rmat(f,N+e,s)*Rmat(e,N,t)
        src,tgt=basis(N),basis(N+e+f); L=N+e+f+4
        for j,lam in enumerate(src):
            Ml=maya(lam,L)
            for i,mu in enumerate(tgt):
                if C[i,j]==0: continue
                if len(Ml ^ maya(mu,L))==2:
                    nz1+=1
                    if sp.simplify(sp.expand(C[i,j].subs(s,1/t)))!=0: surv+=1
    print(f"   (e,f)=({e},{f}): {nz1} nonzero one-bead entries, {surv} SURVIVE on ts=1")
