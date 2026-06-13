#!/usr/bin/env python3
"""
Job B — three-row M_j closed form + sharp J* box.

For lambda=(a,b,c) |- 2m, a>=b>=c>=1, m<=12:
  M_j = <s_lambda, h_1^{2(m-j)} e_2^j>.
  val(j) = j + 2 v2(C(m,j) M_j),  J* = argmin.

(1) Test M_j ?= f^{(a-j, b-j, c-j)}  (3-row dim, rows shortened by j; partition iff c-j>=0).
(2) Independent route: signed S_3 Jacobi-Trudi determinant of trinomials
      M_j = sum_{w in S_3} sgn(w) T_j(lambda + delta - w.delta),  delta=(2,1,0),
      T_j(alpha) = [s^a1 t^a2 u^a3] e_1^{2(m-j)} e_2^j   (trinomial e_2 expansion).
    Cross-checked against the chain engine.
(3) J* census for three-row shapes: |J*|, the sets, the sharp box, first |J*|=4.
"""
from math import factorial
from itertools import permutations
from job_jstar_engine import Mj_chain_all, hook_f, C, v2

# ---------------------------------------------------------------------------
# (2) signed S_3 JT determinant of trinomials -> independent M_j
# ---------------------------------------------------------------------------
def multinomial3(n, a, b, c):
    if a < 0 or b < 0 or c < 0 or a + b + c != n:
        return 0
    return factorial(n) // (factorial(a) * factorial(b) * factorial(c))

def Tj(alpha, j, m):
    """[s^a1 t^a2 u^a3] e_1^{2(m-j)} e_2^j  via trinomial expansion."""
    a1, a2, a3 = alpha
    if min(a1, a2, a3) < 0:
        return 0
    if a1 + a2 + a3 != 2 * m:
        return 0
    A = 2 * (m - j)
    tot = 0
    # e_2^j = sum_{p+q+r=j} (j!/p!q!r!) s^{p+q} t^{p+r} u^{q+r}
    for p in range(j + 1):
        for q in range(j - p + 1):
            r = j - p - q
            coef_e2 = factorial(j) // (factorial(p) * factorial(q) * factorial(r))
            i = a1 - (p + q)
            k = a2 - (p + r)
            l = a3 - (q + r)
            if i < 0 or k < 0 or l < 0 or i + k + l != A:
                continue
            tot += coef_e2 * multinomial3(A, i, k, l)
    return tot

def Mj_determinant(lam, j, m):
    a, b, c = lam
    delta = (2, 1, 0)
    lpd = (a + 2, b + 1, c + 0)
    total = 0
    for w in permutations(range(3)):
        # sign of permutation
        sgn = 1
        for x in range(3):
            for y in range(x + 1, 3):
                if w[x] > w[y]:
                    sgn = -sgn
        wdelta = tuple(delta[i] for i in w)
        alpha = tuple(lpd[t] - wdelta[t] for t in range(3))
        total += sgn * Tj(alpha, j, m)
    return total

# ---------------------------------------------------------------------------
# dimension f^{(a-j,b-j,c-j)}
# ---------------------------------------------------------------------------
def dim_shifted(lam, j):
    shifted = tuple(x - j for x in lam)
    if shifted[-1] < 0:
        return 0
    # weakly decreasing automatically since lam is; drop zeros
    sh = tuple(x for x in shifted if x > 0)
    return hook_f(sh)

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def three_row_partitions(twom):
    """all (a,b,c), a>=b>=c>=1, a+b+c=2m."""
    res = []
    for a in range(1, twom + 1):
        for b in range(1, a + 1):
            c = twom - a - b
            if 1 <= c <= b:
                res.append((a, b, c))
    return res

def analyze3(lam):
    lam = tuple(lam)
    m = sum(lam) // 2
    Ms = Mj_chain_all(lam)  # length m+1
    vals = {}
    for j in range(m + 1):
        Mj = Ms[j]; cj = C(m, j)
        if Mj == 0 or cj == 0:
            continue
        vals[j] = j + 2 * (v2(cj) + v2(Mj))
    Vmin = min(vals.values())
    Jstar = sorted(j for j, vv in vals.items() if vv == Vmin)
    return Ms, vals, Jstar, Vmin

def run():
    # -------- cross-check determinant vs engine, and dim conjecture --------
    print("=" * 78)
    print("PART 1: M_j cross-check (engine vs S3-JT determinant) + dim conjecture")
    print("=" * 78)
    det_fail = 0
    dim_hold = 0
    dim_break = 0
    dim_break_examples = []
    total_pairs = 0
    for m in range(2, 13):
        for lam in three_row_partitions(2 * m):
            Ms, vals, Jstar, Vmin = analyze3(lam)
            for j in range(m + 1):
                Mdet = Mj_determinant(lam, j, m)
                if Mdet != Ms[j]:
                    det_fail += 1
                    if det_fail <= 6:
                        print(f"  DET FAIL lam={lam} j={j}: engine={Ms[j]} det={Mdet}")
                dimv = dim_shifted(lam, j)
                total_pairs += 1
                if Ms[j] == dimv:
                    dim_hold += 1
                else:
                    dim_break += 1
                    if len(dim_break_examples) < 20:
                        dim_break_examples.append((lam, j, Ms[j], dimv))
    print(f"\nTotal (lam,j) pairs (3-row, m<=12): {total_pairs}")
    print(f"Determinant vs engine FAILURES: {det_fail}  (0 = the trinomial determinant IS M_j)")
    print(f"Dimension conjecture M_j = f^(a-j,b-j,c-j):  HOLD {dim_hold}  BREAK {dim_break}")
    if dim_break_examples:
        print("\n  First breaks of M_j = f^(a-j,b-j,c-j):")
        print(f"  {'lam':>12} {'j':>2} {'M_j':>10} {'f^(shift)':>10}")
        for lam, j, mj, dv in dim_break_examples:
            print(f"  {str(lam):>12} {j:>2} {mj:>10} {dv:>10}")

    # characterise WHERE the dim conjecture holds
    print("\n  --- where does M_j = f^(a-j,b-j,c-j) hold? ---")
    # hypotheses: holds for j=0 always (M_0=f^lam); maybe holds while c-j>=0 & ... test
    j0_ok = True
    holds_when_partition = 0; breaks_when_partition = 0
    holds_when_notpart = 0; breaks_when_notpart = 0
    for m in range(2, 13):
        for lam in three_row_partitions(2 * m):
            Ms, *_ = analyze3(lam)
            a, b, c = lam
            for j in range(m + 1):
                dimv = dim_shifted(lam, j)
                is_part = (c - j >= 0)
                if Ms[j] == dimv:
                    if is_part: holds_when_partition += 1
                    else: holds_when_notpart += 1
                else:
                    if is_part: breaks_when_partition += 1
                    else: breaks_when_notpart += 1
            if Ms[0] != hook_f(lam):
                j0_ok = False
    print(f"  M_0 = f^lam always: {j0_ok}")
    print(f"  When c-j>=0 (shift is a partition):  hold {holds_when_partition}  break {breaks_when_partition}")
    print(f"  When c-j<0  (shift not a partition):  hold {holds_when_notpart}  break {breaks_when_notpart}")

    # -------- J* census for three-row shapes --------
    print("\n" + "=" * 78)
    print("PART 2: J* census over three-row lambda, m<=12")
    print("=" * 78)
    from collections import Counter
    size_counter = Counter()
    box_fail = 0
    first_4 = None
    jstar_sets = Counter()
    parity_mixed = 0
    union_offsets = set()
    examples_by_size = {}
    for m in range(2, 13):
        for lam in three_row_partitions(2 * m):
            Ms, vals, Jstar, Vmin = analyze3(lam)
            sz = len(Jstar)
            size_counter[sz] += 1
            off = tuple(j - Jstar[0] for j in Jstar)
            jstar_sets[off] += 1
            union_offsets.update(off)
            # parity
            par = set(j % 2 for j in Jstar)
            if len(par) > 1:
                parity_mixed += 1
            # XOR-closed 2-adic box test on offsets
            offset_set = set(off)
            if not is_xor_box(offset_set):
                box_fail += 1
            if sz == 4 and first_4 is None and len(examples_by_size.get(4, [])) == 0:
                first_4 = (m, lam, Jstar)
            if sz >= 4:
                examples_by_size.setdefault(sz, []).append((m, lam, Jstar))
    print(f"|J*| census (three-row): {dict(sorted(size_counter.items()))}")
    print(f"parity-mixed J* count: {parity_mixed}")
    print(f"XOR-2-adic-box failures on offsets: {box_fail}")
    print(f"union of all J* offsets seen: {sorted(union_offsets)}")
    print(f"distinct J* offset patterns: {dict(sorted(jstar_sets.items()))}")
    if first_4:
        print(f"\nFirst |J*|=4 among three-row shapes: m={first_4[0]} lam={first_4[1]} J*={first_4[2]}")
    for sz in sorted(examples_by_size):
        if sz >= 4:
            ex = examples_by_size[sz][:8]
            print(f"\n|J*|={sz} examples (first {len(ex)}):")
            for m_, lam_, js_ in ex:
                print(f"   m={m_} lam={lam_} J*={js_}")

def is_xor_box(offset_set):
    """Is offset_set an XOR-closed subspace of (Z/2)^inf (as integers), containing 0?"""
    if 0 not in offset_set:
        # offsets are relative to min, so min->0 always present
        return False
    elems = sorted(offset_set)
    # build GF(2) span of nonzero elements, check it equals the set
    basis = []
    for x in elems:
        cur = x
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur != 0:
            basis.append(cur)
            basis.sort(reverse=True)
    # span size
    span = {0}
    for b in basis:
        span |= {s ^ b for s in span}
    return span == offset_set

if __name__ == "__main__":
    run()
