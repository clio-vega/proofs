#!/usr/bin/env python3
"""
Job A (2026-06-16 code session): pin the c=4 three-row boundary 2-adic structure.

lambda=(a,b,4), a+b=2m-4 (a,b SAME parity since c=4 even), P=m-b-4, a=2P+b+4.
Box interior: theta=0, j0=0 (box {0,2,4,6}), need STRICT Delta(b+i)>0 for i=1,2,3,4.

This sweep:
  1. validates the SymPy closed forms M_{b+i} against Murnaghan-Nakayama (mn.Mj);
  2. tabulates v2(N_i) and confirms the memory's CONTENT bounds (0 violations);
  3. computes Delta(b+i) two ways (direct closed form + val(j)-val(0) from MN) and
     records the min slack of Delta(b+i)>0 and the extremal shape;
  4. residue law for v2(M_{b+2}) by (a%8,b%8,P%8): which classes force which v2.

Reuses mn.py from threerow-boundary/.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'threerow-boundary'))
from mn import Mj
from math import comb, gcd
from collections import defaultdict
import sympy as sp

def v2(n):
    n = abs(int(n))
    if n == 0:
        return 99
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k

# ---- closed forms (fitted in c4_forms.py, validated below) ----
A, B = sp.symbols('a b')
M_forms = {
    4: (B-3)*(B+2)*(B+3)*(B+4)/24,
    3: (B-3)*(B+2)*(B+3)*(A*B+4*A-B**2-3*B-8)/24,
    2: (B-3)*(B+2)*(A**2*B**2+7*A**2*B+12*A**2-2*A*B**3-9*A*B**2-21*A*B-20*A
                    + B**4+2*B**3+13*B**2+16*B-8)/48,
    1: (B-3)*(A-B+1)*(A**2*B**3+9*A**2*B**2+26*A**2*B+24*A**2-2*A*B**4-7*A*B**3
                      -13*A*B**2-14*A*B+B**5-2*B**4+17*B**3-4*B**2-84*B-96)/144,
}
# exact integer evaluators (avoid float overflow at large m)
def _int_eval(expr):
    p = sp.Poly(sp.together(expr).as_numer_denom()[0], A, B)
    den = sp.together(expr).as_numer_denom()[1]
    den = int(den)
    terms = [(int(pa), int(pb), int(co)) for (pa, pb), co in p.terms()]
    def f(a, b):
        s = 0
        for pa, pb, co in terms:
            s += co * a**pa * b**pb
        assert s % den == 0, (a, b, s, den)
        return s // den
    return f
M_lam = {i: _int_eval(M_forms[i]) for i in M_forms}
# the N_i deficit factors (the irreducible part)
N_forms = {
    3: A*B+4*A-B**2-3*B-8,
    2: A**2*B**2+7*A**2*B+12*A**2-2*A*B**3-9*A*B**2-21*A*B-20*A+B**4+2*B**3+13*B**2+16*B-8,
    1: A**2*B**3+9*A**2*B**2+26*A**2*B+24*A**2-2*A*B**4-7*A*B**3-13*A*B**2-14*A*B
       +B**5-2*B**4+17*B**3-4*B**2-84*B-96,
}
N_lam = {i: _int_eval(N_forms[i]) for i in N_forms}

c = 4

from math import factorial
def hook_f(lam):
    # hook-length formula for f^lambda (exact, fast)
    lam = [x for x in lam if x > 0]
    rows = len(lam)
    conj = [sum(1 for r in lam if r > col) for col in range(lam[0])]
    n = sum(lam)
    prod = 1
    for i in range(rows):
        for jcol in range(lam[i]):
            hook = (lam[i] - jcol - 1) + (conj[jcol] - i - 1) + 1
            prod *= hook
    return factorial(n) // prod

def box_interior(a, b):
    # c=4 even: a,b same parity. interior gate matching c4_scout: b >= 2c.
    return b >= 2*c

# ============================================================
# 1. Validate closed forms vs MN (small m, exhaustive)
# ============================================================
def validate_forms(mmax=46):
    mism = 0; tested = 0
    for m in range(8, mmax):
        for b in range(c, m):
            a = 2*m - c - b
            if a < b:
                continue
            for i in (1, 2, 3, 4):
                j = b + i
                if j > m:
                    continue
                M_true = Mj((a, b, c), j, m)
                M_form = M_lam[i](a, b)
                tested += 1
                if M_true != M_form:
                    mism += 1
                    if mism <= 8:
                        print(f"   MISMATCH M_{{b+{i}}} (a,b,m)=({a},{b},{m}): MN={M_true} form={M_form}")
    print(f"[1] closed-form vs MN: tested={tested}, mismatches={mism}")
    return mism

# ============================================================
# 2 & 3.  big sweep: content bounds + Delta slacks + residue law
# ============================================================
def big_sweep(mmax=400):
    shapes = 0
    # content-bound violation counters
    viol = defaultdict(int)
    # Delta slack tracking: Delta(b+i) must be > 0 (theta=0). slack = Delta(b+i).
    min_slack = {1: (None, None), 2: (None, None), 3: (None, None), 4: (None, None)}
    delta_dist = {i: defaultdict(int) for i in (1, 2, 3, 4)}
    # residue law for v2(M_{b+2}) and v2(N_i)
    res_Mb2 = defaultdict(set)            # (a%8,b%8,P%8) -> set of v2(M_{b+2})
    v2N = {i: defaultdict(int) for i in (1, 2, 3)}   # v2(N_i) distribution
    v2N_by_bpar = {i: {0: defaultdict(int), 1: defaultdict(int)} for i in (1, 2, 3)}
    mn_delta_mismatch = 0; mn_delta_checked = 0

    for m in range(8, mmax+1):
        for b in range(c, m):
            a = 2*m - c - b
            if a < b:
                continue
            if not box_interior(a, b):
                continue
            P = m - b - c
            if P < 0:
                continue
            shapes += 1

            # val(0); M_0 = f^{(a,b,c)} via hook-length (fast/exact)
            M0 = hook_f((a, b, c))
            val0 = 0 + 2*v2(comb(m, 0)*M0)   # = 2*v2(M0)

            # N_i content checks
            for i in (1, 2, 3):
                Ni = N_lam[i](a, b)
                vN = v2(Ni)
                v2N[i][vN] += 1
                v2N_by_bpar[i][b % 2][vN] += 1
            # memory bounds
            N3 = N_lam[3](a, b); N2 = N_lam[2](a, b); N1 = N_lam[1](a, b)
            if b % 2 == 0 and v2(N3) < 1: viol['N3>=1 (b even)'] += 1
            if v2(N2) < 2: viol['N2>=2'] += 1
            if v2(N1) < 2: viol['N1>=2'] += 1
            if b % 2 == 1 and v2(N1) < 2: viol['N1>=2 (b odd sharp)'] += 1
            if b % 2 == 0 and v2(N1) < 3: viol['N1>=3 (b even)'] += 1

            # Delta(b+i) two ways
            for i in (1, 2, 3, 4):
                j = b + i
                if j > m:
                    continue
                Mform = M_lam[i](a, b)
                if Mform == 0:
                    # M=0 means val=+inf, Delta=+inf, trivially >0
                    continue
                valj_form = j + 2*v2(comb(m, j)*Mform)
                delta = valj_form - val0
                delta_dist[i][delta] += 1
                # track min slack (we need delta>0, i.e. delta>=1; slack=delta)
                cur = min_slack[i][0]
                if cur is None or delta < cur:
                    min_slack[i] = (delta, (a, b, m, P))
                # cross-check against MN val (only small m -- MN is slow)
                if m <= 60:
                    M_mn = Mj((a, b, c), j, m)
                    valj_mn = (j + 2*v2(comb(m, j)*M_mn)) if M_mn != 0 else None
                    mn_delta_checked += 1
                    if M_mn != Mform:
                        mn_delta_mismatch += 1
                # residue law for i=2
                if i == 2:
                    res_Mb2[(a % 8, b % 8, P % 8)].add(v2(Mform))

    print(f"[2/3] big sweep: {shapes} box-interior shapes, m<= {mmax}")
    print(f"   MN Delta cross-check: checked={mn_delta_checked}, mismatches={mn_delta_mismatch}")
    print()
    print("   --- content-bound violations (memory's claims) ---")
    if not viol:
        print("   0 violations of ALL memory content bounds "
              "(N3>=1 b-even, N2>=2, N1>=2, N1>=3 b-even).")
    for k, v in sorted(viol.items()):
        print(f"   VIOLATION {k}: {v}")
    print()
    print("   --- v2(N_i) distributions (overall / b-even / b-odd) ---")
    for i in (1, 2, 3):
        print(f"   N{i}: overall {dict(sorted(v2N[i].items()))}")
        print(f"        b-even {dict(sorted(v2N_by_bpar[i][0].items()))} | "
              f"b-odd {dict(sorted(v2N_by_bpar[i][1].items()))}")
    print()
    print("   --- Delta(b+i) min slack (need Delta>0) and extremal shape ---")
    for i in (1, 2, 3, 4):
        sl, sh = min_slack[i]
        print(f"   Delta(b+{i}): min={sl} at (a,b,m,P)={sh}")
        print(f"        dist (smallest 6): {dict(sorted(delta_dist[i].items())[:6])}")
    print()
    print("   --- v2(M_{b+2}) residue law by (a%8,b%8,P%8) -> set of v2 ---")
    multi = {k: sorted(v) for k, v in res_Mb2.items() if len(v) > 1}
    single = {k: next(iter(v)) for k, v in res_Mb2.items() if len(v) == 1}
    print(f"   classes with a UNIQUE v2(M_b+2): {len(single)} of {len(res_Mb2)}")
    # summarize: which (b%8) forces what; the deficit class
    from collections import Counter
    by_b = defaultdict(set)
    for (a8, b8, P8), vs in res_Mb2.items():
        by_b[b8] |= vs
    print(f"   v2(M_b+2) values seen, grouped by b%8: "
          f"{ {k: sorted(v) for k, v in sorted(by_b.items())} }")
    # worst (largest v2 = most cancellation = closest to deficit) is irrelevant;
    # the *smallest* v2 governs the closest-to-tie. report global min of Delta already above.
    return min_slack

if __name__ == "__main__":
    print("="*70)
    print("Job A: c=4 three-row boundary 2-adic structure")
    print("="*70)
    mm = validate_forms()
    if mm == 0:
        print("   closed forms VALIDATED.\n")
    big_sweep(mmax=400)
