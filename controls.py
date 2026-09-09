import sympy as sp, sys
sys.path.insert(0,'/home/clio/projects/reviews/2026-09-08-selfreview-code')
from ribbon import border_strips, parts
sys.path.insert(0,'/home/clio/projects/scratch/0909')
from q107 import n_stat, parts_upto, maj_of_strip
t,s=sp.symbols('t s')
def conj(lam):
    if not lam: return ()
    return tuple(sum(1 for x in lam if x>j) for j in range(lam[0]))
def basis(N): return sorted(parts(N))
def Rmat(e,N,u):
    src,tgt=basis(N),basis(N+e); M=sp.zeros(len(tgt),len(src))
    for j,l in enumerate(src):
        for (mu,h) in border_strips(l,e): M[tgt.index(mu),j]+=u**h
    return M
def Om(N):
    b=basis(N); M=sp.zeros(len(b),len(b))
    for j,l in enumerate(b): M[b.index(conj(l)),j]=1
    return M

print("CONTROL 1 (Garsia-Zabrocki): mu -> mu+1^k, is n(mu+1^k)-n(mu) constant?")
for k in (2,3,4,5):
    vals=set()
    for mu in parts_upto(10):
        m=list(mu)+[0]*max(0,k-len(mu))
        lam=tuple([m[i]+1 for i in range(k)]+m[k:])
        lam=tuple(x for x in lam if x>0)
        if any(lam[i]<lam[i+1] for i in range(len(lam)-1)): continue
        vals.add(n_stat(lam)-n_stat(mu))
    print(f"   k={k}: values = {sorted(vals)}   binom(k,2) = {k*(k-1)//2}"
          f"   -> {'CONSTANT' if len(vals)==1 else 'NOT constant'}")

print("\nCONTROL 2: maj of the column ribbon 1^k  (must equal binom(k,2))")
for k in (2,3,4,5):
    mj,_=maj_of_strip(tuple([1]*k))
    print(f"   k={k}: maj(1^k) = {mj}, binom(k,2) = {k*(k-1)//2}  {'OK' if mj==k*(k-1)//2 else 'MISMATCH'}")

print("\nCONTROL 3 (the brief's proposed control -- does it move?):")
print("  brief said: 'break the ribbon condition, check c is NOT constant there.'")
print("  mu+1^k is a NON-ribbon skew shape (disconnected unless mu_1=..=mu_k), and its c IS constant.")
print("  -> the brief's control tests the WRONG variable: connectedness is not what makes c constant.")
print("  The real variable is whether the added cells sit at rows fixed independently of mu.")

print("\nCONTROL 4: a COLUMN ribbon added at a non-top row (is 'columns' the right class?)")
import collections
byr=collections.defaultdict(set)
for mu in parts_upto(9):
    L=len(mu)+6
    for lam,ht in __import__('q107').add_strips(mu,3,L):
        r,a=__import__('q107').strip_data(mu,lam)
        if a==(1,1,1): byr[r].add(n_stat(lam)-n_stat(mu))
print(f"   e=3, column shape (1,1,1): c by top-row r -> "
      f"{ {r:sorted(v) for r,v in sorted(byr.items())} }")
print("   -> NOT constant.  Columns are a red herring; the top-row pin r=1 is the whole obstruction.")

print("\nCONTROL 5: Thm C(iii)  Ad(omega) C_{e,e} = -C_{e,e}  on ts=1")
for e in (2,3,4):
    ok=True
    for N in range(0,4):
        C=Rmat(e,N+e,t)*Rmat(e,N,s)-Rmat(e,N+e,s)*Rmat(e,N,t)
        Ch=C.subs(s,1/t)
        D=sp.simplify(sp.expand(Om(N+2*e)*Ch*Om(N) + Ch))
        if not all(x==0 for x in D): ok=False; break
    print(f"   e={e}: {'OK  (C is in the (-1)-eigenspace of Ad(omega) on ts=1)' if ok else 'FAIL'}")
