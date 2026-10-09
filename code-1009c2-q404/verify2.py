"""Fast cross-mechanism verification.  One matrix inverse per n:
   if A has columns = P_mu in the p-basis, then Y = A^{-1} (rows mu, cols rho)."""
import sys, sympy as sp
from itertools import product as iproduct
sys.path.insert(0,'.')
from hl import partitions, HL_P, t
import kf
def flush(*a): print(*a); sys.stdout.flush()
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7

Ymat={}
for n in range(1,NMAX+1):
    parts,Ps=HL_P(n); idx={p:i for i,p in enumerate(parts)}
    A=sp.zeros(len(parts),len(parts))
    for j in range(len(parts)):
        for k,v in Ps[j].items(): A[idx[k],j]=v
    Ainv=A.inv()
    Ymat[n]=(parts,{(parts[i],parts[r]):sp.expand(sp.cancel(Ainv[i,r]))
                    for i in range(len(parts)) for r in range(len(parts))})
    flush('built n=%d (%d partitions)'%(n,len(parts)))

def Yv(mu,rho): return Ymat[sum(mu)][1][(mu,tuple(sorted(rho,reverse=True)))]
def eq(x,y): return sp.expand(x-y)==0

# CHECK1: Y(0) = chi (Murnaghan-Nakayama -- an independent mechanism)
bad=tot=0
for n in range(1,NMAX+1):
    for mu in partitions(n):
        for rho in partitions(n):
            tot+=1
            if Yv(mu,rho).subs(t,0)!=kf.chi(mu,rho): bad+=1
flush('CHECK1  Y^mu_rho(0) = chi^mu_rho (Murnaghan-Nakayama) : %d pairs n<=%d, %d failures'%(tot,NMAX,bad))
badc=0
for n in range(1,NMAX+1):
    for mu in partitions(n):
        for rho in partitions(n):
            if (Yv(mu,rho)+1).subs(t,0)!=kf.chi(mu,rho): badc+=1
flush('  control (Y -> Y+1)                                  : %d of %d fire'%(badc,tot))

# CHECK2: Y(1) = [m_mu] p_rho by direct monomial count
def coeff_m(mu,rho):
    L=len(mu); c=0
    for f in iproduct(range(L),repeat=len(rho)):
        s=[0]*L
        for j,i in enumerate(f): s[i]+=rho[j]
        if tuple(s)==tuple(mu): c+=1
    return c
bad2=tot2=0
for n in range(1,NMAX+1):
    for mu in partitions(n):
        for rho in partitions(n):
            tot2+=1
            if Yv(mu,rho).subs(t,1)!=coeff_m(mu,rho): bad2+=1
flush('CHECK2  Y^mu_rho(1) = [m_mu] p_rho (direct count)      : %d pairs n<=%d, %d failures'%(tot2,NMAX,bad2))

# CHECK3: Theorem C on the two-row x two-part slice
def thmC(a,b,m):
    if m==b: return t**b-t**(b-1)+1+(1 if a==b else 0)
    if 1<=m<=b-1: return (t-1)*t**(b-1)+(t-1)*t**(b-m-1)
    return (t-1)*t**(b-1)
bad3=tot3=0; cells=[]
for n in range(2,NMAX+1):
    for b in range(1,n//2+1):
        a=n-b
        if a<b: continue
        for x in range(1,n//2+1):
            y=n-x; m=min(x,y)
            mu=(a,b); rho=tuple(sorted((x,y),reverse=True))
            tot3+=1; got=Yv(mu,rho); want=thmC(a,b,m)
            if not eq(got,want):
                bad3+=1
                if bad3<6: flush('  FAIL ThmC mu=%s rho=%s got=%s want=%s'%(mu,rho,got,want))
            cells.append((mu,rho,m,sp.expand(got)))
flush('CHECK3  Theorem C closed form                          : %d two-row/two-part cells n<=%d, %d failures'%(tot3,NMAX,bad3))
# control: drop the [a=b] term; predict WHICH cells fire, in advance
pred=[(mu,rho) for mu,rho,m,g in cells if mu[0]==mu[1] and m==mu[1]]
fired=[(mu,rho) for mu,rho,m,g in cells if not eq(g, thmC(mu[0],mu[1],m)-(1 if (mu[0]==mu[1] and m==mu[1]) else 0))]
flush('  control (drop [a=b]): predicted %d cells, fired on %d, same set: %s'%(len(pred),len(fired),sorted(pred)==sorted(fired)))

# CHECK4: the classification of Theorem 6, tested against the ENGINE (not against Theorem C)
sys.path.insert(0,'.')
from exact import unitary_exact
bad4=tot4=0; wrong=[]
for mu,rho,m,g in cells:
    a,b=mu
    pred_unitary = (m!=b) or (b==1) or (b==2 and a>=3)
    ok,_=unitary_exact(g)
    tot4+=1
    if ok!=pred_unitary: bad4+=1; wrong.append((mu,rho,m,g,ok,pred_unitary))
flush('CHECK4  Theorem 6 classification vs exact factorisation: %d cells, %d failures'%(tot4,bad4))
for w in wrong[:5]: flush('   ',w)
nonu=[(mu,rho) for mu,rho,m,g in cells if not unitary_exact(g)[0]]
flush('  non-unitary cells found (n<=%d): %d  -> %s'%(NMAX,len(nonu),nonu))
