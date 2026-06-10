"""
JOB A (arms PROVE) — a closed/2-adic formula for v2(M_j).

M_j = <s_lam, h1^{2(m-j)} e2^j>,  M_0 = f^lam.
val(j) = j + 2 v2(C(m,j) M_j),  J* = argmin val,  mu = min val.

This script:
  (1) tabulates v2(M_j) for ALL lam |- 2m, m<=MMAX, with covariates
      (j-bits, m-bits, 4-core, 4-quotient, f^lam, v2(f), chi_b);
      writes results/v2Mj-table.csv  (one row per (lam,j) with M_j != 0).
  (2) tests candidate closed forms for v2(M_j).
  (3) isolates the step-law decomposition  Dv2(M) = 2^{a-1} - Dv2(binom) on J*.

Run:  python3 jobA_v2Mj.py [MMAX]   (default 12)
"""
import sys, math, csv, os
from collections import Counter, defaultdict
from job1_tie_census import partitions, M_vector, v2, chi_b_vector, analyze
from core_quotient import core_quotient

RESULTS = "/home/clio/projects/code/results"
os.makedirs(RESULTS, exist_ok=True)

def bits(n):
    return bin(n)[2:] if n else "0"

def kummer_carries(a, b):
    """#carries when adding a+b in base 2 = v2(C(a+b, a))."""
    carries = 0; carry = 0
    while a or b or carry:
        s = (a & 1) + (b & 1) + carry
        carry = 1 if s >= 2 else 0
        carries += carry
        a >>= 1; b >>= 1
    return carries

# ---------------------------------------------------------------- table
def build_rows(MMAX):
    rows = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            r = analyze(lam, m)
            Jstar = set(r['Jstar'])
            f = Ms[0]
            v2f = v2(f)
            core, quot = core_quotient(lam, 4)
            quotsizes = tuple(sum(q) for q in quot)
            for j in range(m+1):
                if Ms[j] == 0:
                    continue
                vM = v2(Ms[j]); vb = v2(math.comb(m, j))
                rows.append(dict(
                    m=m, lam="|".join(map(str, lam)), j=j,
                    Mj=Ms[j], v2Mj=vM, v2binom=vb, val=j+2*vb+2*vM,
                    f=f, v2f=(v2f if v2f is not None else -1),
                    jbits=bits(j), mbits=bits(m),
                    in_Jstar=int(j in Jstar),
                    core4="|".join(map(str, core)),
                    quot4=";".join("|".join(map(str, q)) if q else "" for q in quot),
                    quotsize="|".join(map(str, quotsizes)),
                ))
    return rows

# ---------------------------------------------------------------- fits
def fit_tests(MMAX):
    """Test candidate closed forms for v2(M_j) on all shapes m<=MMAX."""
    print(f"\n=== FORMULA FITS for v2(M_j)  (m<={MMAX}) ===")

    # gather (m, lam, j, vM, vf, Ms) tuples
    data = []
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            vf = v2(Ms[0])
            for j in range(m+1):
                if Ms[j] == 0:
                    continue
                data.append((m, lam, j, v2(Ms[j]), vf, Ms))

    N = len(data)
    print(f"data points (M_j != 0): {N}")

    # Candidate 0: monotone floor — is v2(M_j) >= v2(M_0)=v2(f)?
    ge_f = sum(1 for (m,lam,j,vM,vf,Ms) in data if vM >= vf)
    print(f"[C0] v2(M_j) >= v2(f^lam):                 {ge_f}/{N}")

    # Candidate 1: v2(M_j) == v2(f) + (something independent of f)?
    #   test whether  v2(M_j) - v2(f)  depends only on (m, j, 4-quotient)
    #   i.e. is it constant across shapes sharing (m,j,4-core?,4-quotient)?
    by_key = defaultdict(set)
    for (m,lam,j,vM,vf,Ms) in data:
        core, quot = core_quotient(lam, 4)
        key = (m, j, tuple(tuple(q) for q in quot))
        by_key[key].add(vM - vf)
    consistent = sum(1 for v in by_key.values() if len(v) == 1)
    print(f"[C1] v2(M_j)-v2(f) determined by (m,j,4-quot): "
          f"{consistent}/{len(by_key)} keys single-valued")

    # Candidate 2: Kummer — does v2(M_j) match carries of some (j + g) base 2?
    #   For each shape, find if there is a SHIFT g (0<=g) with
    #   v2(M_j) = kummer_carries(j, g)  for all j.  Report fraction of shapes.
    shapes_fit = 0; shapes_tot = 0
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            vMs = {j: v2(Ms[j]) for j in range(m+1) if Ms[j] != 0}
            shapes_tot += 1
            found = False
            for g in range(0, 4*m+1):
                if all(vMs[j] == kummer_carries(j, g) for j in vMs):
                    found = True; break
            if found:
                shapes_fit += 1
    print(f"[C2] per-shape Kummer-shift v2(M_j)=carries(j,g): {shapes_fit}/{shapes_tot} shapes")

    # Candidate 3: Lucas product over base-2 digits of j against m.
    #   v2 of a Lucas product C(m,j) mod 2 is 0/positive; not directly v2(M).
    #   Instead test: is v2(M_j) additive over the SET BITS of j?
    #   i.e. v2(M_j) == sum_{a: bit a of j set} delta(m,lam,a) for some per-bit deltas.
    shapes_add = 0; shapes_tot2 = 0
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms = M_vector(lam, m)
            vMs = {j: v2(Ms[j]) for j in range(m+1) if Ms[j] != 0}
            if 0 not in vMs:
                continue
            shapes_tot2 += 1
            base = vMs[0]
            # per-bit delta from single-bit j = 2^a present in domain
            nbits = m.bit_length()
            delta = {}
            ok = True
            for a in range(nbits+1):
                ja = 1 << a
                if ja in vMs:
                    delta[a] = vMs[ja] - base
            # predict all j
            for j, vM in vMs.items():
                pred = base + sum(delta.get(a, None) if (j>>a)&1 else 0
                                  for a in range(nbits+1))
                if any((j>>a)&1 and a not in delta for a in range(nbits+1)):
                    ok = False; break
                if pred != vM:
                    ok = False; break
            if ok:
                shapes_add += 1
    print(f"[C3] v2(M_j) additive over set-bits of j (base v2(M_0)): "
          f"{shapes_add}/{shapes_tot2} shapes")

    return data

# ---------------------------------------------------------------- step law
def steplaw_decomp(MMAX):
    """On every involution pair {j, j+2^a} of every tie's J*, isolate the
    M-part and binomial-part of the step law:
        Dval = (j+2^a + 2 v2(C M)_{j+2^a}) - (j + 2 v2(C M)_j) = 0  on J*.
    =>  2^a + 2[ Dv2(binom) + Dv2(M) ] = 0  where D = (j+2^a)-part minus j-part.
    =>  Dv2(M) = -2^{a-1} - Dv2(binom).
    Quantify how often the binomial part is 0 (M carries everything) vs not.
    """
    from job2_mechanism import val_pieces, jstar_of, box_generators
    print(f"\n=== STEP-LAW DECOMPOSITION on J* pairs  (m<={MMAX}) ===")
    pairs=0; ok=0
    bin_zero=0; bin_nonzero=0
    Mcarries_all=0
    dv2M_counter=Counter(); dv2bin_counter=Counter()
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms,pj,pbin,pM,valv = val_pieces(lam,m)
            J,mu = jstar_of(valv)
            if len(J) < 2:
                continue
            j0,S = box_generators(J)
            if not S:
                continue
            Jset=set(J)
            for a in S:
                g=1<<a
                for j in J:
                    if ((j-j0)>>a)&1:
                        continue
                    jp=j+g
                    if jp not in Jset:
                        continue
                    pairs+=1
                    dv2bin = v2(math.comb(m,jp)) - v2(math.comb(m,j))
                    dv2M   = v2(Ms[jp]) - v2(Ms[j])
                    # step law: dv2M = -2^{a-1} - dv2bin   (since (jp-j)=2^a)
                    if dv2M == -(2**(a-1)) - dv2bin:
                        ok+=1
                    dv2M_counter[(a,dv2M)]+=1
                    dv2bin_counter[(a,dv2bin)]+=1
                    if dv2bin==0:
                        bin_zero+=1
                        if dv2M == -(2**(a-1)):
                            Mcarries_all+=1
                    else:
                        bin_nonzero+=1
    print(f"involution pairs: {pairs}")
    print(f"  step law dv2M = -2^(a-1) - dv2bin holds: {ok}/{pairs}")
    print(f"  pairs with binomial part ZERO (M carries all): {bin_zero}/{pairs}")
    print(f"    of those, dv2M == -2^(a-1) exactly:          {Mcarries_all}/{bin_zero}")
    print(f"  pairs with binomial part NONZERO:              {bin_nonzero}/{pairs}")
    print(f"  dv2M distribution by generator a: ", dict(sorted(dv2M_counter.items())))
    print(f"  dv2bin distribution by generator a:", dict(sorted(dv2bin_counter.items())))

# ---------------------------------------------------------------- S bound
def S_bound(MMAX):
    from job2_mechanism import val_pieces, jstar_of, box_generators
    print(f"\n=== S subset {{1,2,3}} bound  (m<={MMAX}) ===")
    maxS_by_m = defaultdict(int)
    maxgen_overall = 0
    viol = []
    maxS_vs_size=Counter()
    for m in range(1, MMAX+1):
        for lam in partitions(2*m):
            Ms,pj,pbin,pM,valv = val_pieces(lam,m)
            J,mu = jstar_of(valv)
            if len(J) < 2:
                continue
            j0,S = box_generators(J)
            if not S:
                continue
            mx=max(S)
            maxS_by_m[m]=max(maxS_by_m[m], mx)
            maxgen_overall=max(maxgen_overall, mx)
            maxS_vs_size[(len(J), mx)]+=1
            if mx > 3:
                viol.append((m,lam,J,S))
    print(f"  max generator exponent overall: {maxgen_overall}")
    print(f"  max(S) by m:", dict(sorted(maxS_by_m.items())))
    print(f"  (|J*|, max generator) distribution:", dict(sorted(maxS_vs_size.items())))
    if viol:
        print(f"  *** {len(viol)} VIOLATIONS of S subset {{1,2,3}} (generator >= 16):")
        for v in viol[:20]:
            print("    ", v)
    else:
        print(f"  NO generator >= 16: S subset {{1,2,3}} CONFIRMED to m={MMAX}.")

if __name__ == "__main__":
    MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    rows = build_rows(MMAX)
    with open(f"{RESULTS}/v2Mj-table.csv","w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"wrote results/v2Mj-table.csv ({len(rows)} rows, m<={MMAX})")
    fit_tests(MMAX)
    steplaw_decomp(MMAX)
    S_bound(MMAX)
