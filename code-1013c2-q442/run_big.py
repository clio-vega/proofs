"""Larger confirmation + NEGATIVE CONTROLS.  Engine B (Laurent) only where cheap."""
import sys, time
sys.path.insert(0,'/home/clio/projects/proofs/code-1013c2-q442')
from poly import *
from q442 import *

def eq(p,q): return p==q
def nonneg(p): return all(c>0 for c in p.values())

DOMAINS = [(1,0,n) for n in range(2,9)] + [(1,1,n) for n in range(2,9)] + \
          [(1,2,n) for n in range(3,9)] + [(1,3,n) for n in range(4,9)] + \
          [(2,0,n) for n in range(2,7)] + [(2,1,n) for n in range(2,7)] + \
          [(2,2,n) for n in range(3,7)] + [(3,0,n) for n in range(2,5)] + \
          [(3,1,n) for n in range(2,5)]
RUN_B = {(k,l,n) for (k,l,n) in DOMAINS if (k==1 and n<=6) or (k==2 and n<=4) or (k==3 and n<=3)}

print("DOMAIN LIST (%d triples): %s" % (len(DOMAINS), DOMAINS))
print("engine B (Warnaar Laurent) run on %d of them: %s" % (len(RUN_B), sorted(RUN_B)))
print()
hdr=(f"{'(k,ell,n)':>12} {'w':>3} | {'A=B':>5} | {'MAIN A=Ainf':>11} | {'scr [y0]=[k0]':>13} "
     f"{'diff>0':>7} {'diff=ow+Bi':>11} {'n>=w':>5} | {'NC1 A=Ctot':>10} {'NC2 A=Abd':>9} "
     f"{'NC3 A=k0':>8} | {'#tdeg':>5} {'t1gate':>6}")
print(hdr); print("-"*len(hdr))
S=dict(ab=0,abN=0,main=0,scr_fire=0,scr_blind=0,pos=0,pred=0,nc1=0,nc2=0,nc3=0,t1=0,tot=0,
       scr_pred=0)
t0=time.time()
for (k,ell,n) in DOMAINS:
    w=2*ell+2
    A=engineA(k,n); C=engineC_sums(k,n,w)
    dpU=engineC(k,n,w,bounded=False)
    Ctot={}
    for lam,p in dpU.items(): Ctot=add(Ctot,p)         # NC1: no sign rule at all
    k0=C['kappa0']; diff=sub(A,k0); pred=add(C['overwide'],C['Binf'])
    ab=None
    if (k,ell,n) in RUN_B:
        B,G1,td=engineB(k,ell,n); ab=eq(A,B); S['abN']+=1; S['ab']+=ab
        S['t1']+=eq(G1,sub(C['A'],C['B'])); ntd=len(td)
    else:
        ntd=None
    S['tot']+=1
    main=eq(A,C['Ainf']); S['main']+=main
    fire=not eq(A,k0)
    S['scr_fire']+=fire; S['scr_blind']+= (not fire)
    S['scr_pred']+= (fire == (n>=w))
    if fire: S['pos']+= (1 if nonneg(diff) else 0)
    S['pred']+=eq(diff,pred)
    nc1=eq(A,Ctot); nc2=eq(A,C['A']); nc3=eq(A,k0)
    S['nc1']+= (not nc1); S['nc2']+= (not nc2); S['nc3']+= (not nc3)
    print(f"{str((k,ell,n)):>12} {w:>3} | {str(ab):>5} | {str(main):>11} | {str(not fire):>13} "
          f"{str(nonneg(diff) if diff else '-'):>7} {str(eq(diff,pred)):>11} {str(n>=w):>5} | "
          f"{str(nc1):>10} {str(nc2):>9} {str(nc3):>8} | {str(ntd):>5} ")
print("-"*len(hdr))
T=S['tot']
print(f"\nDOMAIN: {T} triples (k,ell,n), k in 1..3, ell in 0..3, n in 2..8 (k=1), 2..6 (k=2), 2..4 (k=3).")
print(f"  engine A == engine B (Warnaar [t^0], Laurent)        : {S['ab']}/{S['abN']}")
print(f"  t=1 gate G(1)(x..)^k == A_alpha - B_alpha            : {S['t1']}/{S['abN']}")
print(f"  MAIN  det_k(f_{{i-j}}-f_{{i+j}}) == sum_a A^inf_a x^a  : {S['main']}/{T}")
print(f"  SCREEN fires ([y=0] != [kappa=0])                    : {S['scr_fire']}/{T}  "
      f"(blind: {S['scr_blind']})")
print(f"  screen fires  <=>  n >= w                            : {S['scr_pred']}/{T}")
print(f"  difference strictly positive where nonzero           : {S['pos']}/{S['scr_fire']}")
print(f"  difference == overwide-even + B^inf  (structure)      : {S['pred']}/{T}")
print(f"\nNEGATIVE CONTROLS (each must FAIL, i.e. count = {T}/{T} 'differs'):")
print(f"  NC1  A != (all lambda, no sign rule)                 : {S['nc1']}/{T}")
print(f"  NC2  A != A_alpha (bounded/cylindric count)          : {S['nc2']}/{T}  "
      f"(expected < {T}: equal where the bound is inactive)")
print(f"  NC3  A != [kappa=0]                                  : {S['nc3']}/{T}  "
      f"(= the screen; expected < {T}: equal where n<w)")
print(f"\nelapsed {time.time()-t0:.1f}s")
