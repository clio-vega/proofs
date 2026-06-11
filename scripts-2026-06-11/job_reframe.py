"""
Reframing test for even-|J*|.

G = A + (-pi_bar) B,   A = sum_t C(m,2t) M_{2t} (-2i)^t,  B = sum_t C(m,2t+1) M_{2t+1} (-2i)^t.
Parity-even J*  ->  governed by A = a(-2i),  a(s)=sum_t alpha_t s^t,  alpha_t = C(m,2t) M_{2t}.
  v_pi(alpha_t (-2i)^t) = 2(v2(alpha_t)+t).  So |J*| = #{t : v2(alpha_t)+t = min} = on-line pts of
  the HONEST 2-adic Newton polygon of a(s) on the slope-(-1) line.
Parity-odd similarly with b(s)=sum_t C(m,2t+1) M_{2t+1} s^t.

Verify: this count == |J*| from job_jstar.  Also dump the box generators (differences) for even cases.
"""
import sys
import symfunc as sf
import job_jstar as JJ
from collections import Counter

def v2(n):
    return JJ.v2(n)

def half_count(coeffs):
    # coeffs[t] = alpha_t (integers, some maybe 0).  return #minimizers of v2(alpha_t)+t and the minimizer t-set
    best=None; args=[]
    for t,a in enumerate(coeffs):
        if a==0: continue
        val=v2(a)+t
        if best is None or val<best:
            best=val; args=[t]
        elif val==best:
            args.append(t)
    return len(args), args, best

def analyze(lam,m):
    d=JJ.analyze(lam,m)
    M=d['M']
    from math import comb
    # even half: alpha_t = C(m,2t) M_{2t}
    alpha=[comb(m,2*t)*M[2*t] for t in range((m)//2+1)]
    # odd half: beta_t = C(m,2t+1) M_{2t+1}
    beta=[comb(m,2*t+1)*M[2*t+1] for t in range((m-1)//2+1)]
    parity = d['V']%2
    if parity==0:
        cnt,args,best=half_count(alpha)
    else:
        cnt,args,best=half_count(beta)
    return d, cnt, args, parity

def main():
    mmax=int(sys.argv[1]) if len(sys.argv)>1 else 7
    mismatch=[]
    gens=Counter()
    for m in range(1,mmax+1):
        for lam in sf.partitions(2*m):
            d,cnt,args,parity=analyze(lam,m)
            if cnt!=d['absJ']:
                mismatch.append((lam,m,d['absJ'],cnt,parity))
            if d['absJ']>=2:
                Js=d['Jstar']; diffs=tuple(Js[k]-Js[0] for k in range(len(Js)))
                gens[diffs]+=1
    print(f"reframing mismatches (half-poly 2-adic count != |J*|): {len(mismatch)}")
    for x in mismatch[:10]: print("   MM:",x)
    print(f"J* offset-patterns (Js - min): {dict(sorted(gens.items(),key=lambda kv:(len(kv[0]),kv[0])))}")

if __name__=="__main__":
    main()
