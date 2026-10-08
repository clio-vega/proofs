from fractions import Fraction
from math import gcd, factorial
from lattice import *
from entrywise import aut_group, act, canon
import sys
def lcm(a,b): return a*b//gcd(a,b) if a and b else 0
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 8

def analyse(lam):
    k=len(lam); dl=d_of(lam); G=aut_group(lam)
    SP=[canon(pi) for pi in set_partitions(k)]
    fib={}
    for pi in SP: fib.setdefault(coarsen(lam,pi),[]).append(pi)
    out={}
    for nu,F in fib.items():
        N=sum(mu_hat0(pi) for pi in F); ent=Fraction(N,dl)
        seen=set(); orbs=[]
        for pi in F:
            if pi in seen: continue
            orb={act(w,pi) for w in G}; seen|=orb
            orbs.append((dl//len(orb), mu_hat0(pi)))
        bound=1
        for stab,m in orbs: bound=lcm(bound, stab//gcd(stab,abs(m)) if m else 1)
        out[nu]=(ent,bound,orbs)
    return dl,out

print("=== (A) is the ROW maximum attained only at the diagonal?  (rows with d_lam>1 only) ===")
print("    rows with d_lam=1 are vacuous: every denominator is 1 = d_lam there.")
rows_v=0; rows_f=0; only_diag=0; more=0; examples=[]
for n in range(1,NMAX+1):
    for lam in partitions(n):
        dl,out=analyse(lam)
        if dl==1: rows_v+=1; continue
        rows_f+=1
        att=[nu for nu,(e,b,o) in out.items() if e.denominator==dl]
        assert lam in att or tuple(lam) in att, ('diagonal must attain', lam, att)
        if len(att)==1: only_diag+=1
        else:
            more+=1
            if len(examples)<6: examples.append((n,lam,dl,sorted(att)))
print("    vacuous rows (d_lam=1)        : %d"%rows_v)
print("    FALSIFIABLE rows (d_lam>1)    : %d"%rows_f)
print("    attained ONLY at the diagonal : %d"%only_diag)
print("    attained elsewhere too        : %d"%more)
for e in examples: print("      n=%d lam=%s d=%d attaining nu: %s"%e)

print()
print("=== (B) first column: entry at nu=(n) is (-1)^{k-1}(k-1)!/d_lam, denominator d_lam/gcd(d_lam,(k-1)!) ===")
bad=0; chk=0; nontriv=0
for n in range(1,NMAX+1):
    for lam in partitions(n):
        dl,out=analyse(lam); k=len(lam)
        ent,bound,orbs=out[(n,)]
        pred=Fraction((-1)**(k-1)*factorial(k-1), dl)
        predden=dl//gcd(dl,factorial(k-1))
        chk+=1
        if ent!=pred or ent.denominator!=predden: bad+=1
        if dl>1 and predden!=dl: nontriv+=1   # the cases where the first column DROPS below d_lam
print("    rows checked                            : %d"%chk)
print("    rows where the drop is strict (d>1, den<d): %d   <-- falsifiable: the formula predicts a PROPER divisor"%nontriv)
print("    mismatches                              : %d"%bad)

print()
print("=== (C) PLANTED CONTROL on the orbit-stabiliser bound ===")
print("    The bound is 'den divides lcm_O(|Stab_O|/gcd(|Stab_O|,mu_O))'.  Plant the")
print("    off-by-duality error |Stab_O| := |O| (orbit size instead of stabiliser order).")
print("    If the test cannot tell these apart it is reading nothing.")
for label,use_orbit in (("TRUE stabiliser",False),("PLANTED: orbit size",True)):
    viol=0; fals=0
    for n in range(1,NMAX+1):
        for lam in partitions(n):
            k=len(lam); dl=d_of(lam); G=aut_group(lam)
            SP=[canon(pi) for pi in set_partitions(k)]
            fib={}
            for pi in SP: fib.setdefault(coarsen(lam,pi),[]).append(pi)
            for nu,F in fib.items():
                N=sum(mu_hat0(pi) for pi in F); den=Fraction(N,dl).denominator
                seen=set(); bound=1; norb=0
                for pi in F:
                    if pi in seen: continue
                    orb={act(w,pi) for w in G}; seen|=orb; norb+=1
                    q = len(orb) if use_orbit else dl//len(orb)
                    m=mu_hat0(pi)
                    bound=lcm(bound, q//gcd(q,abs(m)) if m else 1)
                if norb>1: fals+=1
                if den and bound%den: viol+=1
    print("      %-22s : %d violations  (out of %d support entries, %d multi-orbit)"%(label,viol,396 if NMAX==8 else -1,fals))
