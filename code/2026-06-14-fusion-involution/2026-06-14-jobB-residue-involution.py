#!/usr/bin/env python3
"""
JOB B (06-14 system cycle) -- confirm the even-|J*| fpf involution is the
TOP-GENERATOR toggle, pairing by leading-unit residue mod pi^2, NOT Sawin adjacency.

SETUP (verified identity, see job_jstar_engine + census script):
  G_lambda(i) = <s_lambda,(p2+pi e2)^m> = Sum_j C(m,j) M_j (i-1)^j,   pi=1+i,
  v_pi(term_j) = val(j) = j + 2 v2(C(m,j) M_j).
  J* = argmin_j val(j).  Leading unit  u_j = C(m,j) M_j (i-1)^j / pi^V   (Gaussian integer, v_pi=0).
  Even-|J*| <=> the leading layer C = Sum_{j in J*} u_j cancels (C = |J*| mod pi).

RESIDUE MOD pi^2:  pi^2 = (1+i)^2 = 2i, and the ideal (2i)=(2), so
  Z[i]/(pi^2) = Z[i]/(2)  has 4 elements; residue(u) = (Re(u) mod 2) + (Im(u) mod 2) i.
  A unit (v_pi=0) has odd norm => exactly one of Re,Im odd => residue in {1, i}.

DISCRIMINATOR for a pairing {j, j'}:  sum u_j + u_j' should have v_pi >= 2 (cancel to next layer).
  - residue-equal partners (r, r):  u_j+u_j' = 2r + (higher) ,  v_pi(2 r) = 2  => CANCELS.
  - residue-opposite (1, i):  u_j+u_j' ~ 1+i = pi,  v_pi = 1            => does NOT cancel.

CLAIM (to confirm):  the fpf involution is XOR-with-the-TOP-generator of the J* box
  (for the |J*|=4 family J*={3,5,7,9}=3+{0,2,4,6} this is j <-> j+4), and it pairs
  residue-EQUAL partners on every tie, while Sawin's adjacent j<->j+2 pairs residue-opposite.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # bundled job_jstar_engine
sys.path.insert(0, '/home/clio/projects/code')
from sympy import I, expand, Integer, re, im
from math import comb
from collections import Counter, defaultdict
from job_jstar_engine import Mj_chain_all, C, v2, analyze, partitions

# ----------------------------------------------------------------------
def vpi(z):
    """pi=(1+i)-adic valuation of nonzero Gaussian integer z; v_pi(z)=v2(norm(z))."""
    z = expand(z)
    a = int(re(z)); b = int(im(z))
    if a == 0 and b == 0:
        return None
    return v2(a*a + b*b)

def residue_pi2(u):
    """Residue of Gaussian integer u mod pi^2 = mod 2, as a label in {'1','i','0','1+i'}."""
    u = expand(u)
    a = int(re(u)) % 2; b = int(im(u)) % 2
    return {(0,0):'0', (1,0):'1', (0,1):'i', (1,1):'1+i'}[(a, b)]

def leading_units(lam):
    """Return (units dict {j:u_j}, V, J*, m). u_j = C(m,j)M_j(i-1)^j / pi^V, a Gaussian integer."""
    A = analyze(lam)
    m, J, V, M = A['m'], A['Jstar'], A['V'], A['M']
    piV = (1 + I) ** V
    units = {}
    for j in J:
        term = Integer(comb(m, j)) * M[j] * (I - 1) ** j
        u = expand(term / piV)           # already a Gaussian integer (u_j = +-odd or +-odd*i)
        a, b = re(u), im(u)
        assert a == int(a) and b == int(b), f"non-integer unit at lam={lam}, j={j}: {u}"
        units[j] = Integer(int(a)) + Integer(int(b)) * I
    return units, V, J, m

def xor_top_generator(J):
    """The XOR-with-top-generator involution on the affine 2-adic box J (offsets from min).
    Returns dict j -> sigma(j). Top generator g = largest 2^a appearing in the offset GF(2) span."""
    j0 = J[0]
    offs = [j - j0 for j in J]
    # GF(2) span generators of the offsets (offsets are multiples of 2 here, but treat generally)
    basis = []
    for x in sorted(offs):
        cur = x
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur); basis.sort(reverse=True)
    if not basis:
        return {j: j for j in J}
    g = max(basis)                      # TOP generator (largest power of two)
    sigma = {}
    for j in J:
        sigma[j] = j0 + ((j - j0) ^ g)  # XOR offset with top generator
    return sigma

# ----------------------------------------------------------------------
def residue_table(lam):
    units, V, J, m = leading_units(lam)
    res = {j: residue_pi2(units[j]) for j in J}
    sigma = xor_top_generator(J)
    return units, V, J, m, res, sigma

def check_pairing_cancels(units, pairs):
    """For a list of (j,j') pairs, return list of (pair, sum, v_pi(sum), cancels?)."""
    out = []
    for (j, jp) in pairs:
        s = expand(units[j] + units[jp])
        v = vpi(s) if s != 0 else 'oo'
        cancels = (s == 0) or (isinstance(v, int) and v >= 2)
        out.append(((j, jp), s, v, cancels))
    return out

def matching_from_involution(sigma):
    """Turn an fpf involution dict into a set of unordered pairs; None if it has a fixed point."""
    pairs = []
    seen = set()
    for j in sigma:
        if sigma[j] == j:
            return None
        if j in seen:
            continue
        pairs.append((j, sigma[j])); seen.add(j); seen.add(sigma[j])
    return pairs

# ----------------------------------------------------------------------
def run():
    print("="*86)
    print("JOB B: even-|J*| fpf involution = XOR-with-top-generator; residue mod pi^2 discriminator")
    print("="*86)

    # ---- (1) the |J*|=4 residue tables: c=3 family (a,b,3) with a odd, plus extras ----
    print("\n[1] |J*|=4 residue tables (the discriminating case: 3 possible pairings, 1 cancels)\n")
    # gather some |J*|=4 shapes (three-row c=3 'both generators' family) + minimal witness
    quad_shapes = []
    for a in range(3, 30):
        for b in range(3, a+1):
            lam = (a, b, 3)
            if sum(lam) % 2: continue
            if sum(lam)//2 > 16: continue
            A = analyze(lam)
            if A['absJstar'] == 4:
                quad_shapes.append(lam)
    # also look for |J*|=4 in c=2 / non-three-row up to m<=14 by scanning partitions
    for n in range(4, 29, 2):
        for lam in partitions(n):
            if len(lam) < 2: continue
            A = analyze(lam)
            if A['absJstar'] == 4 and lam not in quad_shapes:
                quad_shapes.append(lam)
    quad_shapes = quad_shapes[:8]

    for lam in quad_shapes:
        units, V, J, m, res, sigma = residue_table(lam)
        print(f"  lam={lam}  m={m}  V={V}  J*={J}")
        for j in J:
            print(f"     j={j}: u_j={units[j]!s:>14}   residue mod pi^2 = {res[j]}")
        # three pairings of a 4-set
        a0,a1,a2,a3 = J
        pairings = {
            "TOP-gen (XOR top)": matching_from_involution(sigma),
            "adjacent (Sawin j<->j+2)": [(a0,a1),(a2,a3)],
            "outer/cross": [(a0,a3),(a1,a2)],
        }
        for name, pairs in pairings.items():
            if pairs is None:
                print(f"     {name:28s}: involution has a FIXED POINT (not fpf)"); continue
            chk = check_pairing_cancels(units, pairs)
            allc = all(c for *_, c in chk)
            detail = "; ".join(f"{{{p[0]},{p[1]}}}:v_pi(sum)={v}{'OK' if c else 'x'}"
                                for (p, s, v, c) in chk)
            tag = "<-- CANCELS (the involution)" if allc and name.startswith("TOP") else ("ALL CANCEL" if allc else "")
            print(f"     {name:28s}: {detail}   {tag}")
        print()

    # ---- (2) fpf on EVERY tie + the residue-mod-4 law that drives the pairing ----
    print("="*86)
    print("[2] XOR-top involution fpf on EVERY tie; residue law (j=j' mod4 <=> residue-equal)")
    print("="*86)
    print("  NOTE: leading-layer cancellation (evenness) needs only mod pi: every unit = 1 = i mod pi")
    print("        (since 1-i = -i*pi), so ANY pairing sums to 0 mod pi once |J*| is even.")
    print("        The mod pi^2 refinement (v_pi(sum)>=2) is a |J*|>=4 phenomenon; for |J*|=2 the")
    print("        single pair reaches mod pi^2 iff its generator is 4 (residue-equal), not 2.\n")
    # scan: three-row families up to m<=16 (cheap) + all partitions for small m
    tie_shapes = []
    for a in range(2, 32):
        for b in range(1, a+1):
            for c in range(0, b+1):
                lam = tuple(x for x in (a,b,c) if x>0)
                if sum(lam) % 2 or sum(lam)//2 > 16 or sum(lam) < 2: continue
                A = analyze(lam)
                if A['absJstar'] % 2 == 0 and A['absJstar'] >= 2:
                    tie_shapes.append(lam)
    for n in range(2, 23, 2):       # all partitions of n<=22 (m<=11)
        for lam in partitions(n):
            A = analyze(lam)
            if A['absJstar'] % 2 == 0 and A['absJstar'] >= 2 and lam not in tie_shapes:
                tie_shapes.append(lam)
    tie_shapes = sorted(set(tie_shapes), key=lambda L:(sum(L), L))

    fpf_fail = []
    res_law_fail = []          # residue(u_j)==residue(u_j') XOR (j==j' mod 4): should be all-equal
    modpi_fail = []            # leading layer C=sum units must be 0 mod pi (v_pi(C)>=1)
    modpi2_quad_fail = []      # for |J*|>=4: top-gen pairs must cancel mod pi^2
    size_ctr = Counter()
    for lam in tie_shapes:
        A = analyze(lam)
        size_ctr[A['absJstar']] += 1
        units, V, J, m = leading_units(lam)
        sigma = xor_top_generator(J)
        pairs = matching_from_involution(sigma)
        if pairs is None:
            fpf_fail.append((lam, J)); continue
        res = {j: residue_pi2(units[j]) for j in J}
        # residue law: equal residue  <=>  j == j' (mod 4)
        for jx in J:
            for jy in J:
                if (res[jx] == res[jy]) != ((jx - jy) % 4 == 0):
                    res_law_fail.append((lam, jx, jy, res[jx], res[jy]))
        # leading-layer (mod pi) cancellation = evenness
        Csum = expand(sum(units.values()))
        if not (Csum == 0 or (vpi(Csum) is not None and vpi(Csum) >= 1)):
            modpi_fail.append((lam, J, vpi(Csum)))
        # mod pi^2 per-pair for |J*|>=4
        if A['absJstar'] >= 4:
            chk = check_pairing_cancels(units, pairs)
            if not all(c for *_, c in chk):
                modpi2_quad_fail.append((lam, J))

    print(f"  even-|J*| tie shapes checked: {len(tie_shapes)}   |J*| distribution: {dict(sorted(size_ctr.items()))}")
    print(f"  (a) XOR-top fixed-point failures (not fpf):                    {len(fpf_fail)}")
    print(f"  (b) residue law  [res equal <=> j=j' mod4]  failures:          {len(res_law_fail)}")
    print(f"  (c) leading-layer mod-pi cancellation (evenness) failures:     {len(modpi_fail)}")
    print(f"  (d) |J*|>=4: top-gen pairs cancel mod pi^2 -- failures:        {len(modpi2_quad_fail)}")
    for f in fpf_fail[:5]: print("      fpf:", f)
    for f in res_law_fail[:5]: print("      reslaw:", f)
    for f in modpi_fail[:5]: print("      modpi:", f)
    for f in modpi2_quad_fail[:5]: print("      modpi2:", f)
    verdict = not (fpf_fail or res_law_fail or modpi_fail or modpi2_quad_fail)
    print(f"\n  VERDICT [2]: involution fpf, residue-mod-4 law holds, evenness cancels mod pi,")
    print(f"               and on every |J*|>=4 tie the TOP-generator pair cancels mod pi^2: "
          f"{'CONFIRMED' if verdict else 'FALSIFIED'}")

    # ---- (3) Sawin adjacency comparison on |J*|=4 ----
    print("\n" + "="*86)
    print("[3] Sawin adjacent j<->j+2 vs top-gen j<->j+4 on |J*|=4 ties")
    print("="*86)
    saw_ok = 0; saw_fail = 0; top_ok = 0
    for lam in quad_shapes:
        units, V, J, m = leading_units(lam)
        a0,a1,a2,a3 = J
        sigma = xor_top_generator(J)
        top_pairs = matching_from_involution(sigma)
        adj_pairs = [(a0,a1),(a2,a3)]
        top_c = all(c for *_, c in check_pairing_cancels(units, top_pairs))
        adj_c = all(c for *_, c in check_pairing_cancels(units, adj_pairs))
        top_ok += top_c
        if adj_c: saw_ok += 1
        else: saw_fail += 1
        print(f"  lam={lam}: top-gen cancels={top_c}, Sawin-adjacent cancels={adj_c}")
    print(f"\n  top-gen cancels on {top_ok}/{len(quad_shapes)};  Sawin-adjacent cancels on {saw_ok}/{len(quad_shapes)}")
    print(f"  => Sawin adjacency {'FALSIFIED' if saw_fail else 'matches'} ; top-generator toggle is the involution.")

if __name__ == "__main__":
    run()
