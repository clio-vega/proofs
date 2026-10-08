from fractions import Fraction
from math import gcd
from lattice import *
from entrywise import aut_group, act, canon
import sys
def lcm(a,b): return a*b//gcd(a,b) if a and b else 0
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 8
print("ENTRYWISE BOUND, with the falsifiable set counted FIRST")
print("  (tightness is an identity when the fibre is a single orbit: entry = mu_O/|Stab_O|.")
print("   The claim can only fail when the fibre carries >= 2 orbits.)")
print()
hdr="  %-4s %-7s %-7s %-9s %-9s %-7s %-7s"%("n","support","1-orbit","MULTI","viol","tight","loose")
print(hdr); print("  "+"-"*62)
tot={'support':0,'single':0,'multi':0,'viol':0,'tight':0,'loose':0}
loose_list=[]
for n in range(1,NMAX+1):
    c={'support':0,'single':0,'multi':0,'viol':0,'tight':0,'loose':0}
    for lam in partitions(n):
        k=len(lam); dl=d_of(lam); G=aut_group(lam)
        SP=[canon(pi) for pi in set_partitions(k)]
        fib={}
        for pi in SP: fib.setdefault(coarsen(lam,pi),[]).append(pi)
        for nu,F in fib.items():
            N=sum(mu_hat0(pi) for pi in F); ent=Fraction(N,dl); den=ent.denominator
            seen=set(); bound=1; norb=0
            for pi in F:
                if pi in seen: continue
                orb={act(w,pi) for w in G}; seen|=orb; norb+=1
                stab=dl//len(orb); m=mu_hat0(pi)
                bound=lcm(bound, stab//gcd(stab,abs(m)) if m else 1)
            c['support']+=1
            if norb==1: c['single']+=1
            else: c['multi']+=1
            if den and bound%den: c['viol']+=1
            if den==bound: c['tight']+=1
            else:
                c['loose']+=1
                if norb>1: loose_list.append((n,lam,nu,norb,den,bound,ent))
    print("  %-4d %-7d %-7d %-9d %-9d %-7d %-7d"%(n,c['support'],c['single'],c['multi'],c['viol'],c['tight'],c['loose']))
    for k2 in tot: tot[k2]+=c[k2]
print("  "+"-"*62)
print("  %-4s %-7d %-7d %-9d %-9d %-7d %-7d"%("all",tot['support'],tot['single'],tot['multi'],tot['viol'],tot['tight'],tot['loose']))
print()
print("  VERDICT on the falsifiable set (multi-orbit fibres), n <= %d:"%NMAX)
print("    falsifiable instances : %d"%tot['multi'])
print("    bound violated        : %d"%tot['viol'])
print("    bound loose (strict)  : %d"%len(loose_list))
print("    bound tight           : %d"%(tot['multi']-len(loose_list)))
print()
print("  the loose instances (cross-orbit cancellation), all of them:")
for e in loose_list: print("    n=%d lam=%s nu=%s orbits=%d exact_den=%d bound=%d entry=%s"%e)
