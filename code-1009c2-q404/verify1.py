"""Cross-mechanism verification of the Green-polynomial engine, and of Theorem C."""
import sys, sympy as sp
from itertools import product as iproduct
sys.path.insert(0,'.')
from hl import partitions, Y_def, t
import kf

NMAX = int(sys.argv[1]) if len(sys.argv)>1 else 6

def flush(*a):
    print(*a); sys.stdout.flush()

# ---- CHECK 1: Y(0) = chi^lam_rho  (Murnaghan-Nakayama: a different mechanism) ----
bad=0; tot=0
for n in range(1,NMAX+1):
    for lam in partitions(n):
        for rho in partitions(n):
            tot+=1
            if sp.simplify(Y_def(lam,rho).subs(t,0) - kf.chi(lam,rho))!=0:
                bad+=1
                if bad<4: flush('  FAIL t=0',lam,rho)
flush('CHECK1  Y^lam_rho(0) = chi^lam_rho   : %d (lam,rho) pairs n<=%d, %d failures'%(tot,NMAX,bad))

# ---- CHECK 2: Y(1) = [m_lam] p_rho (direct monomial count: another mechanism) ----
def coeff_m(lam,rho):
    L=len(lam); cnt=0
    for f in iproduct(range(L),repeat=len(rho)):
        s=[0]*L
        for j,i in enumerate(f): s[i]+=rho[j]
        if tuple(s)==tuple(lam): cnt+=1
    return cnt
bad2=0; tot2=0
for n in range(1,NMAX+1):
    for lam in partitions(n):
        for rho in partitions(n):
            tot2+=1
            if sp.simplify(Y_def(lam,rho).subs(t,1)-coeff_m(lam,rho))!=0:
                bad2+=1
                if bad2<4: flush('  FAIL t=1',lam,rho,Y_def(lam,rho).subs(t,1),coeff_m(lam,rho))
flush('CHECK2  Y^lam_rho(1) = [m_lam]p_rho  : %d (lam,rho) pairs n<=%d, %d failures'%(tot2,NMAX,bad2))

# ---- CHECK 2b: positive control for CHECK2 -- perturb Y and confirm the test fires ----
bad2c=0; tot2c=0
for n in range(1,NMAX+1):
    for lam in partitions(n):
        for rho in partitions(n):
            tot2c+=1
            if sp.simplify((Y_def(lam,rho)+1).subs(t,1)-coeff_m(lam,rho))!=0: bad2c+=1
flush('CHECK2 control (Y -> Y+1)             : %d pairs, %d failures  [must equal %d]'%(tot2c,bad2c,tot2c))

# ---- CHECK 3: Theorem C closed form on the two-row slice ----
def thmC(a,b,x,y):
    m=min(x,y)
    if m==b:   return t**b - t**(b-1) + 1 + (1 if a==b else 0)
    if 1<=m<=b-1: return (t-1)*t**(b-1) + (t-1)*t**(b-m-1)
    return (t-1)*t**(b-1)
bad3=0; tot3=0; rows=[]
for n in range(2,NMAX+1):
    for b in range(1,n//2+1):
        a=n-b
        if a<b: continue
        for x in range(1,n//2+1):
            y=n-x
            lam=(a,b) if b>0 else (a,)
            lam=tuple(p for p in (a,b) if p>0)
            rho=tuple(sorted((x,y),reverse=True))
            Yv=sp.expand(Y_def(lam,rho)); Cv=sp.expand(thmC(a,b,x,y))
            tot3+=1
            ok = sp.simplify(Yv-Cv)==0
            if not ok:
                bad3+=1
                if bad3<6: flush('  FAIL ThmC lam=%s rho=%s  Y=%s  C=%s'%(lam,rho,Yv,Cv))
            rows.append((lam,rho,min(x,y),Yv,ok))
flush('CHECK3  Theorem C closed form         : %d (lam,rho) on two-row x two-part slice n<=%d, %d failures'%(tot3,NMAX,bad3))

# ---- CHECK 3 control: perturb Theorem C's a==b indicator and confirm it fires ----
def thmC_bad(a,b,x,y):
    m=min(x,y)
    if m==b:   return t**b - t**(b-1) + 1          # DROP the [a=b] term
    if 1<=m<=b-1: return (t-1)*t**(b-1) + (t-1)*t**(b-m-1)
    return (t-1)*t**(b-1)
predicted = sum(1 for (lam,rho,m,Yv,ok) in rows if len(lam)==2 and lam[0]==lam[1] and m==lam[1])
fired=0
for n in range(2,NMAX+1):
    for b in range(1,n//2+1):
        a=n-b
        if a<b: continue
        for x in range(1,n//2+1):
            y=n-x
            lam=tuple(p for p in (a,b) if p>0); rho=tuple(sorted((x,y),reverse=True))
            if sp.simplify(Y_def(lam,rho)-thmC_bad(a,b,x,y))!=0: fired+=1
flush('CHECK3 control (drop [a=b])           : fired on %d pairs, predicted %d (the a=b, min(rho)=b cells)'%(fired,predicted))
