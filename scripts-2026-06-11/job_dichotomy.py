"""
Verify the dichotomy:  |J*| odd  =>  v_pi(G) = V exactly (G != 0).
Compute v_pi(G) from the exact Gaussian integer G (norm), compare to V and |J*| parity.
Also measure pruning power: fraction of shapes with |J*| odd.
"""
import sys
from sympy import Integer, I, expand
import symfunc as sf
from collections import Counter
import job_jstar as JJ

def vpi_int_norm(G):
    # v_pi(G) = v2(N(G)), N = a^2+b^2
    G = complex(G)
    a=int(round(G.real)); b=int(round(G.imag))
    if a==0 and b==0: return None
    N=a*a+b*b
    return JJ.v2(N)

def Gexact(lam,m):
    slam=sf.schur_p(lam)
    psi=sf.add(sf.power(sf.h1(),2), {k:v*(I*(1+I)) for k,v in sf.e_n(2).items()})
    val=expand(sf.inner(slam, sf.power(psi,m)))
    return val

def main():
    mmax=int(sys.argv[1]) if len(sys.argv)>1 else 7
    tot=0; oddJ=0
    bad_odd=[]      # |J*| odd but v_pi(G) != V
    even_vanish=Counter()
    for m in range(1,mmax+1):
        for lam in sf.partitions(2*m):
            tot+=1
            d=JJ.analyze(lam,m)
            G=Gexact(lam,m)
            vG=vpi_int_norm(G)
            if d['absJ']%2==1:
                oddJ+=1
                # claim: vG == V  (and G != 0)
                if vG is None or vG!=d['V']:
                    bad_odd.append((lam,m,d['absJ'],d['V'],vG))
            else:
                # even |J*|: record whether vanishes and excess vG-V
                exc = 'ZERO' if vG is None else vG-d['V']
                even_vanish[exc]+=1
    print(f"total shapes m<={mmax}: {tot}")
    print(f"|J*| odd: {oddJ}  ({100*oddJ/tot:.1f}% pruned by dichotomy)")
    print(f"odd-|J*| with v_pi(G) != V  (would refute lemma): {len(bad_odd)}")
    for x in bad_odd[:10]: print("   BAD:",x)
    print(f"even-|J*| excess (v_pi(G)-V) distribution: {dict(sorted(even_vanish.items(),key=str))}")

if __name__=="__main__":
    main()
