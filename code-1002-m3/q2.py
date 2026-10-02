"""Q2: is beta = sum_{nu in slice} beta_{w(nu)} log-concave in real cylindric instances?
(beta LC ==> G LC by the tail lemma, so this is a STRENGTHENING of condition (A).)
Report census counts.  Also report (A) itself, computed WITHOUT the radial shortcut."""
import gen, layers, sys
from collections import Counter

def run(m, nmax, dmax, report_witness=6):
    st=Counter(); wit_beta=[]; wit_A=[]
    for n in range(m, nmax+1):
        for (mu,lam) in gen.pairs(n,m,dmax):
            d=sum(lam)-sum(mu)
            for b in range(0,d+1):
                off,co = gen.slice_sum(mu,lam,n,m,b)
                if not co: continue
                st['slices']+=1
                # (A): direct PF2 test on the slice sum
                if gen.is_pf2(co): st['A_ok']+=1
                else:
                    st['A_fail']+=1
                    if len(wit_A)<report_witness: wit_A.append((n,mu,lam,b,off,co))
                # beta: radial decrement of G.  Needs G symmetric -- assert it.
                if co != co[::-1]:
                    st['G_not_symmetric']+=1
                    continue
                bet = layers.beta_of(layers.radial(co))
                if min(bet) < 0: st['beta_negative']+=1; continue
                if layers.is_lc(bet): st['beta_lc']+=1
                else:
                    st['beta_notlc']+=1
                    if len(wit_beta)<report_witness: wit_beta.append((n,mu,lam,b,co,bet))
                if layers.is_concave_pospart(bet): st['beta_concave_pospart']+=1
    return st, wit_beta, wit_A

for (m,nmax,dmax) in [(2,10,12),(3,10,12),(4,9,11),(5,9,10)]:
    st,wb,wa = run(m,nmax,dmax)
    print(f"m={m}, n<={nmax}, d<={dmax}: {st['slices']} nonzero slices")
    print(f"   (A) PF2:        ok {st['A_ok']}  FAIL {st['A_fail']}")
    print(f"   G symmetric:    non-symmetric {st['G_not_symmetric']}")
    print(f"   beta >= 0:      negative {st['beta_negative']}")
    print(f"   beta LC:        ok {st['beta_lc']}  FAIL {st['beta_notlc']}")
    print(f"   beta = c_+ (concave): {st['beta_concave_pospart']}")
    for x in wa: print("   (A) WITNESS", x)
    for x in wb: print("   beta WITNESS", x)
    sys.stdout.flush()
