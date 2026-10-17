"""Close the derivation loop:
 (1) reproduce the PAPER's measurement -- engine's type-B object (WB/E^C) collapses with rho^C.
     This is what the paper's 48/48 actually measured.  If it reproduces, then combined with
     the measured E^B/E^C ratio it PREDICTS rho^B = (1-q)^2/((1-q w_a)(1-q w_b)), which
     eb_rho.py then confirmed independently.  Two routes to the same multiplier.
 (2) the true type-B multiplier at N=4 (wider domain; |M| takes three distinct values).
 PLANTS (same block as the theorem -- the multiplier): P1 rho->rho(1+q), P2 rho->1."""
import sympy as sp, itertools
import q484_engine as E
from q484_engine import q, wvars
from eb_audit import Atilde_B_true

def restrictH(f,w,a,b): return sp.cancel(sp.together(f.subs(w[b-1],1/w[a-1])))

def run(N, which, label):
    w=wvars(N); tot=nd=0
    hit={'claim':0,'P1':0,'P2':0}
    for a,b in itertools.combinations(range(1,N+1),2):
        wa=w[a-1]
        rhoC = sp.cancel((1-q)**2/((1-q*wa**2)*(1-q*wa**-2)))
        rhoB = sp.cancel((1-q)**2/((1-q*wa)*(1-q/wa)))
        rho  = rhoC if which=='engineB' else rhoB
        for M in E.subsets(N):
            if not (a in M and b in M): continue
            Mp=tuple(sorted(set(M)-{a,b}))
            for I in E.subsets(N):
                if which=='engineB':
                    AM=E.Atilde(N,w,M,I,X='B'); AMp=E.Atilde(N,w,Mp,I,X='B')
                else:
                    AM=Atilde_B_true(N,w,M,I); AMp=Atilde_B_true(N,w,Mp,I)
                AH,ApH=restrictH(AM,w,a,b),restrictH(AMp,w,a,b)
                tot+=1
                if AH==0 and ApH==0: continue
                nd+=1
                for tag,rr in (('claim',rho),('P1',sp.cancel(rho*(1+q))),('P2',sp.Integer(1))):
                    if sp.simplify(sp.cancel(AH-rr*ApH))==0: hit[tag]+=1
    print(f"  {label} N={N}: pairs {tot}, non-degenerate {nd};  claim {hit['claim']}/{nd}, "
          f"P1 {hit['P1']}/{nd} (must be 0), P2 {hit['P2']}/{nd} (must be 0)")

print("=== (1) reproduce the paper: engine type-B object (WB/E^C) vs rho^C ===")
for N in (2,3): run(N,'engineB','engine-B vs rho^C')
print("=== (2) TRUE type-B (WB/E^B) vs rho^B, N=4 ===")
run(4,'trueB','true-B vs rho^B')
