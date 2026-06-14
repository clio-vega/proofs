#!/usr/bin/env python3
"""
Job B items 2 & 3:
  (2) Census the c=3 J* box: generator set S subset {2,4}, congruence rule in (a,b) mod 2^k.
  (3) Sawin adjacent-pair involution test on the leading-pi coefficient at (9,6,3) and other |J*|=4.

KEY IDENTITY (derived, verified below):
  G_lambda(i) = <s_lambda,(p_2+pi e_2)^m> = Sum_j C(m,j) M_j (i-1)^j,   pi = 1+i,
  because p_2 = h_1^2 - 2 e_2 so p_2 + pi e_2 = h_1^2 + (i-1) e_2, and
  <s_lambda,(h_1^2 + x e_2)^m> = Sum_j C(m,j) M_j x^j.
  v_pi(term_j) = j + 2 v_2(C(m,j) M_j) = val(j)   [since i-1 ~ pi (v_pi=1), v_pi(2)=2].
  J* = argmin_j val(j);  G=0 requires the leading-pi layer C = Sum_{j in J*} unit_j to cancel,
  which needs |J*| even (C = |J*| mod pi).
"""
from sympy import symbols, I, expand, simplify, Integer, nsimplify, re, im
from math import comb
from collections import Counter, defaultdict
from job_jstar_engine import Mj_chain_all, C, v2

def analyze(lam):
    lam=tuple(x for x in lam if x>0); m=sum(lam)//2
    M=Mj_chain_all(lam); vals={}
    for j in range(m+1):
        if M[j]==0 or C(m,j)==0: continue
        vals[j]=j+2*(v2(C(m,j))+v2(M[j]))
    V=min(vals.values()); J=sorted(j for j,v in vals.items() if v==V)
    return M,vals,J,V,m

# ---------- verify the key identity G_lambda(i)=Sum_j C(m,j)M_j(i-1)^j and its v_pi ----------
def Glam_i_via_Mj(lam):
    M,vals,J,V,m=analyze(lam)
    s=Integer(0)
    for j in range(m+1):
        s+= comb(m,j)*M[j]*(I-1)**j
    return expand(s), V

def vpi(z):
    """pi=(1+i)-adic valuation of a Gaussian integer z (z!=0). v_pi(z)=v_2(N(z)) where N=norm."""
    z=expand(z)
    a=int(re(z)); b=int(im(z))
    if a==0 and b==0: return None
    norm=a*a+b*b
    return v2(norm)   # v_pi(z) = v_2(|z|^2)

def leading_units(lam):
    """Return list of (j, unit_j) for j in J*, unit_j = C(m,j)M_j(i-1)^j / pi^V (a Gaussian integer)."""
    M,vals,J,V,m=analyze(lam)
    piV=(1+I)**V
    units=[]
    for j in J:
        term=Integer(comb(m,j))*M[j]*(I-1)**j
        u=expand(term/piV)
        u=nsimplify(u)
        units.append((j, u))
    return units, V, J, m

def run():
    print("="*78); print("VERIFY identity + v_pi(G)=min val unless cancellation"); print("="*78)
    for lam in [(4,3,3),(9,6,3),(6,3,3),(5,4,3),(8,5,3),(2,2),(3,1)]:
        if sum(lam)%2: continue
        G,V=Glam_i_via_Mj(lam)
        vp = vpi(G) if G!=0 else 'oo(G=0)'
        M,vals,J,Vv,m=analyze(lam)
        print(f"  lam={lam}: G_lam(i)={G}, v_pi(G)={vp}, min val(J*)={V}, J*={J}, |J*|={len(J)}")
    print("  -> v_pi(G) > min val  iff leading layer cancels iff |J*| even (necessary).")

    # ---------- (2) c=3 J* box census ----------
    print("\n"+"="*78); print("(2) c=3 J* BOX CENSUS: (a,b,3), generator set S, offsets"); print("="*78)
    rows=[]
    size_ctr=Counter(); gen_ctr=Counter()
    by_S=defaultdict(list)
    for a in range(3, 34):
        for b in range(3, a+1):
            c=3
            if (a+b+c)%2!=0: continue
            lam=(a,b,c); m=(a+b+c)//2
            if m>18: continue
            M,vals,J,V,m=analyze(lam)
            off=tuple(j-J[0] for j in J)   # offsets from min
            # generators of XOR box: offsets should be {0} + subspace of {2,4,6}=2*{1,2,3}
            half=tuple(o//2 for o in off)  # in {0,1,2,3}
            S=xor_generators(set(half))
            size_ctr[len(J)]+=1; gen_ctr[tuple(sorted(S))]+=1
            by_S[tuple(sorted(S))].append((a,b,J[0],J))
            rows.append((a,b,m,J,off,tuple(sorted(S))))
    print(f"  |J*| sizes: {dict(sorted(size_ctr.items()))}")
    print(f"  generator-sets S (in half-units, *2 for actual offset): {dict(gen_ctr)}")
    # first |J*|=4 (S={1,2} i.e. gens 2&4)
    print("\n  Shapes with BOTH generators (|J*|=4, S={1,2} -> offsets {0,2,4,6}):")
    for (a,b,j0,J) in by_S.get((1,2),[])[:12]:
        print(f"    (a,b,3)=({a},{b},3): J*={J}")
    # generator 4 alone (S={2}) -> offsets {0,4}
    print("\n  Shapes with generator 4 only (S={2} -> offsets {0,4}):")
    for (a,b,j0,J) in by_S.get((2,),[])[:12]:
        print(f"    ({a},{b},3): J*={J}")

    # congruence rule: tabulate (a mod 2^k, b mod 2^k) -> S for k=2,3
    print("\n  Congruence hunt: does (a,b) mod 2^k determine S?")
    for k in [2,3,4]:
        mod=2**k
        table=defaultdict(set)
        for (a,b,m,J,off,S) in rows:
            table[(a%mod,b%mod)].add(S)
        well_defined=all(len(v)==1 for v in table.values())
        print(f"    mod {mod}: S determined by (a,b) mod {mod}? {well_defined}"
              f"  ({sum(1 for v in table.values() if len(v)>1)} ambiguous classes)")
        if well_defined:
            print(f"      RULE (a%{mod},b%{mod}) -> S:")
            for key in sorted(table):
                print(f"        {key}: {list(table[key])[0]}")

    # ---------- (3) Sawin adjacent-pair involution test ----------
    print("\n"+"="*78); print("(3) SAWIN adjacent-pair involution on leading-pi units"); print("="*78)
    targets=[(9,6,3)] + [ (a,b,3) for (a,b,j0,J) in by_S.get((1,2),[]) ][:5]
    seen=set()
    for lam in targets:
        if lam in seen: continue
        seen.add(lam)
        units,V,J,m=leading_units(lam)
        print(f"\n  lam={lam}, m={m}, J*={J}, V={V}")
        for j,u in units:
            print(f"     j={j}: unit = C(m,j)M_j(i-1)^j/pi^V = {u}")
        # Sawin: sort J*, pair consecutive (j0,j1),(j2,j3),...; check each pair sums to 0 (cancels)
        uvals=[u for (j,u) in units]
        Js=[j for (j,u) in units]
        print("     -- consecutive-pairing cancellation (Sawin adjacent): --")
        all_cancel=True
        for p in range(0,len(uvals)-1,2):
            ssum=expand(uvals[p]+uvals[p+1])
            cancels=(ssum==0)
            all_cancel=all_cancel and cancels
            print(f"        pair (j={Js[p]},j={Js[p+1]}): {uvals[p]} + {uvals[p+1]} = {ssum}  cancel={cancels}")
        # also test total sum and mod-pi
        tot=expand(sum(uvals))
        print(f"     total leading sum C = {tot}  (=0 means full leading-layer cancellation)")
        print(f"     consecutive-pairs all cancel? {all_cancel}")

def xor_generators(half_set):
    """Given a set of integers (half-offsets) containing 0, return GF(2) basis (generators)."""
    elems=sorted(half_set)
    basis=[]
    for x in elems:
        cur=x
        for b in basis:
            cur=min(cur,cur^b)
        if cur!=0:
            basis.append(cur); basis.sort(reverse=True)
    return set(basis)

if __name__=="__main__":
    run()
