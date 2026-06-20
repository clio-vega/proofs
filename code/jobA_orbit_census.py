#!/usr/bin/env python3
"""Job A (2026-06-20 CODE) -- ORBIT-SIZE CENSUS of the 2-adic toggle, OWED TO RICK.

For three-row interior shapes (a,b,4) and (a,b,5) we proved the argmin set
    J* = argmin_j  val(j),   val(j) = j + 2 v2( C(m,j) M_j ),
        M_j = <s_lambda, h1^{2(m-j)} e2^j>   (Murnaghan-Nakayama, abacus chi engine)
is a SINGLE-generator affine 2-adic box (|S|<=1), so |J*| in {1,2} per shape.

Rick (emails UID 411, 412) wants the ORBIT-SIZE CENSUS: over an interior range,
how many shapes give |J*|=1 (rigid) vs |J*|=2 (toggle), broken down by (a,b) mod 4,
and the TIE LOCUS as a fraction (his "17/17" analogue): of the shapes where the
mod-4 gate makes a tie POSSIBLE, what fraction actually tie -- all-or-nothing per
residue class, or fractional?

GROUND TRUTH = MN (mn.py). Every |J*| is ALSO recomputed from the proved closed
form and the two are asserted equal (0 mismatches reported).

  c=4 closed form: M_j = C(2(m-j),b-j)(a-b+1)Q4 / [24 (a+5-j) prod_{i=1..5... no}1..4(b+i-j)]
  c=5 closed form: M_j = C(2(m-j),b-j)(a-b+1)Q5 / [120(a+6-j) prod_{i=1..5}(b+i-j)]

Predicted tie gates (from the interior proofs):
  c=4: tie {0,2}  <=>  (a,b) = (2,2) or (1,3) (mod 4)         [j0 = 0 always]
  c=5: a even -> j0=0, tie {0,2} <=> (a,b)=(0,1) (mod 4)
       a odd  -> j0=3, tie {3,5} <=> (a,b)=(3,0) (mod 4)
"""
import sys, time, pickle
sys.path.insert(0, '/home/clio/projects/code/threerow-c3')
from mn import Mj
from math import comb
import sympy as sp

def v2(n):
    n = abs(int(n))
    if n == 0:
        return float('inf')
    return (n & -n).bit_length() - 1

# ---------------------------------------------------------------------------
# Closed forms (independently MN-validated below before use).
# ---------------------------------------------------------------------------
a_s, b_s, j_s = sp.symbols('a b j', integer=True)

# Q4 -- octic, lifted from the c=4 interior note (jobA_c4_interior_verify.py).
Q4 = (a_s**4*b_s**4 + 6*a_s**4*b_s**3 - a_s**4*b_s**2 - 54*a_s**4*b_s - 72*a_s**4
 + 10*a_s**3*b_s**4 - 4*a_s**3*b_s**3*j_s**2 - 8*a_s**3*b_s**3*j_s + 60*a_s**3*b_s**3
 - 24*a_s**3*b_s**2*j_s - 10*a_s**3*b_s**2 + 28*a_s**3*b_s*j_s**2 + 80*a_s**3*b_s*j_s
 - 540*a_s**3*b_s + 24*a_s**3*j_s**2 + 192*a_s**3*j_s - 720*a_s**3 + 23*a_s**2*b_s**4
 - 12*a_s**2*b_s**3*j_s**2 - 48*a_s**2*b_s**3*j_s + 138*a_s**2*b_s**3 + 6*a_s**2*b_s**2*j_s**4
 - 12*a_s**2*b_s**2*j_s**3 + 30*a_s**2*b_s**2*j_s**2 - 144*a_s**2*b_s**2*j_s - 23*a_s**2*b_s**2
 - 18*a_s**2*b_s*j_s**4 + 60*a_s**2*b_s*j_s**3 - 6*a_s**2*b_s*j_s**2 + 504*a_s**2*b_s*j_s
 - 1242*a_s**2*b_s - 72*a_s**2*j_s**3 + 72*a_s**2*j_s**2 + 1080*a_s**2*j_s - 1656*a_s**2
 - 34*a_s*b_s**4 + 16*a_s*b_s**3*j_s**2 + 8*a_s*b_s**3*j_s - 204*a_s*b_s**3
 - 6*a_s*b_s**2*j_s**4 + 36*a_s*b_s**2*j_s**3 - 30*a_s*b_s**2*j_s**2 + 48*a_s*b_s**2*j_s
 + 34*a_s*b_s**2 - 4*a_s*b_s*j_s**6 + 48*a_s*b_s*j_s**5 - 226*a_s*b_s*j_s**4 + 492*a_s*b_s*j_s**3
 - 638*a_s*b_s*j_s**2 + 112*a_s*b_s*j_s + 1836*a_s*b_s + 12*a_s*j_s**6 - 144*a_s*j_s**5
 + 732*a_s*j_s**4 - 1800*a_s*j_s**3 + 1752*a_s*j_s**2 - 984*a_s*j_s + 2448*a_s - 120*b_s**4
 + 48*b_s**3*j_s**2 + 240*b_s**3*j_s - 720*b_s**3 - 12*b_s**2*j_s**4 - 24*b_s**2*j_s**3
 - 60*b_s**2*j_s**2 + 672*b_s**2*j_s + 120*b_s**2 + 8*b_s*j_s**6 - 96*b_s*j_s**5 + 524*b_s*j_s**4
 - 1224*b_s*j_s**3 + 1076*b_s*j_s**2 - 2880*b_s*j_s + 6480*b_s + j_s**8 - 28*j_s**7 + 298*j_s**6
 - 1672*j_s**5 + 5305*j_s**4 - 9244*j_s**3 + 9084*j_s**2 - 8928*j_s + 8640)

_d5 = pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl', 'rb'))
Q5 = sp.sympify(_d5['Q5sym'])
# sympify makes fresh symbols named a,b,j -- remap them onto our integer symbols
Q5 = Q5.subs({s: {'a': a_s, 'b': b_s, 'j': j_s}[s.name] for s in Q5.free_symbols})

def _int_eval(poly):
    terms = [(tuple(int(e) for e in mono), int(co))
             for mono, co in sp.Poly(sp.expand(poly), a_s, b_s, j_s).terms()]
    def f(av, bv, jv):
        return sum(co * av**ea * bv**eb * jv**ej for (ea, eb, ej), co in terms)
    return f

Q4_eval = _int_eval(Q4)
Q5_eval = _int_eval(Q5)

def M_closed(c, av, bv, jv, m):
    """Closed form M_j for c in {4,5} as an exact integer."""
    if bv - jv < 0:
        # prefactor C(2(m-j), b-j) vanishes
        return 0
    N = 2 * (m - jv)
    if c == 4:
        num = comb(N, bv - jv) * (av - bv + 1) * Q4_eval(av, bv, jv)
        den = 24 * (av + 5 - jv)
        for i in range(1, 5):
            den *= (bv + i - jv)
    else:  # c == 5
        num = comb(N, bv - jv) * (av - bv + 1) * Q5_eval(av, bv, jv)
        den = 120 * (av + 6 - jv)
        for i in range(1, 6):
            den *= (bv + i - jv)
    q, r = divmod(num, den)
    assert r == 0, f"non-integer closed form c={c} ({av},{bv}) j={jv}"
    return q

# ---------------------------------------------------------------------------
def jstar(c, av, bv, m, use_closed):
    """Return sorted argmin set J* over interior j = 0..b."""
    best = None
    Js = []
    for jv in range(0, bv + 1):
        M = M_closed(c, av, bv, jv, m) if use_closed else Mj((av, bv, c), jv, m)
        if M == 0:
            continue
        val = jv + 2 * v2(comb(m, jv) * M)
        if best is None or val < best:
            best = val
            Js = [jv]
        elif val == best:
            Js.append(jv)
    return Js, best

# Predicted tie gates --------------------------------------------------------
def tie_possible(c, av, bv):
    """Does the mod-4 gate permit a tie for this shape?"""
    ar, br = av % 4, bv % 4
    if c == 4:
        return (ar, br) in {(2, 2), (1, 3)}
    else:  # c == 5
        if av % 2 == 0:
            return (ar, br) == (0, 1)
        else:
            return (ar, br) == (3, 0)

# ---------------------------------------------------------------------------
def run_census(c, MAXM):
    print("=" * 72)
    print(f"CENSUS  c={c},  interior  a>b>{c},  m<= {MAXM}")
    print("=" * 72)
    t0 = time.time()
    # tables
    by_class = {}      # (a%4,b%4) -> [count|J*|=1, count|J*|=2]
    j0_by_class = {}   # (a%4,b%4) -> set of observed j0 = min J*
    mism = 0
    n_shapes = 0
    n_tie = 0
    n_tie_possible = 0
    n_tie_in_possible = 0
    bad_box = 0        # J* not within predicted single-gen box
    examples_tie = {}
    for m in range(c + 2, MAXM + 1):
        n = 2 * m
        for av in range(c + 1, n):
            bv = n - c - av
            if bv <= c or bv >= av:        # interior strict a > b > c
                continue
            n_shapes += 1
            J_mn, V = jstar(c, av, bv, m, use_closed=False)
            J_cf, _ = jstar(c, av, bv, m, use_closed=True)
            if J_mn != J_cf:
                mism += 1
                if mism <= 8:
                    print(f"  MISMATCH MN vs closed c={c} ({av},{bv},{c}) m={m}: "
                          f"MN J*={J_mn} closed J*={J_cf}")
            J = J_mn
            cls = (av % 4, bv % 4)
            slot = by_class.setdefault(cls, [0, 0])
            slot[0 if len(J) == 1 else 1] += 1
            j0_by_class.setdefault(cls, set()).add(min(J))
            if len(J) >= 2:
                n_tie += 1
                if cls not in examples_tie:
                    examples_tie[cls] = (av, bv, m, J)
            poss = tie_possible(c, av, bv)
            if poss:
                n_tie_possible += 1
                if len(J) >= 2:
                    n_tie_in_possible += 1
            # box sanity: |J*|<=2 and J* an arithmetic toggle j0,j0+2
            if len(J) > 2 or (len(J) == 2 and J[1] - J[0] != 2):
                bad_box += 1
                print(f"  BOX VIOLATION c={c} ({av},{bv},{c}) m={m}: J*={J}")
    # report ------------------------------------------------------------
    print(f"\n  shapes swept = {n_shapes}   (elapsed {time.time()-t0:.1f}s)")
    print(f"  MN vs closed-form |J*| MISMATCHES = {mism}")
    print(f"  single-generator box violations (|J*|>2 or step!=2) = {bad_box}")
    print(f"  total toggles |J*|=2 = {n_tie}   rigid |J*|=1 = {n_shapes - n_tie}")
    print()
    print("  orbit-size census by (a%4, b%4):")
    print("    (a%4,b%4) | #|J*|=1 | #|J*|=2 |  j0 seen | tie-gate?")
    for cls in sorted(by_class):
        r1, r2 = by_class[cls]
        gate = "TIE" if (cls in {(2,2),(1,3)} and c==4) else \
               ("TIE" if (c==5 and ((cls==(0,1)) or (cls==(3,0)))) else "")
        print(f"      {cls}   |  {r1:5d}  |  {r2:5d}  |  {sorted(j0_by_class[cls])} | {gate}")
    print()
    # tie locus fraction (Rick's 17/17 analogue)
    if n_tie_possible:
        print(f"  TIE LOCUS (Rick's abundance analogue):")
        print(f"    shapes where mod-4 gate permits a tie : {n_tie_possible}")
        print(f"    of those that ACTUALLY tie            : {n_tie_in_possible}")
        frac = n_tie_in_possible / n_tie_possible
        print(f"    fraction = {n_tie_in_possible}/{n_tie_possible} = {frac:.4f}")
        print(f"    all-or-nothing per residue class?  "
              f"{'YES (=1.0)' if n_tie_in_possible==n_tie_possible else 'NO (fractional)'}")
    # consistency: every tie must lie in a gate class
    stray = n_tie - n_tie_in_possible
    print(f"    ties OUTSIDE the predicted gate classes = {stray}  (expect 0)")
    if examples_tie:
        print("  minimal tie witnesses per class:")
        for cls in sorted(examples_tie):
            print(f"      {cls}: (a,b,m)={examples_tie[cls][:3]} J*={examples_tie[cls][3]}")
    return {
        'c': c, 'n_shapes': n_shapes, 'n_tie': n_tie, 'mism': mism,
        'bad_box': bad_box, 'by_class': by_class,
        'n_tie_possible': n_tie_possible, 'n_tie_in_possible': n_tie_in_possible,
    }


if __name__ == "__main__":
    MAXM = 60
    # First validate the closed forms against MN on a cheap range (belt & braces)
    print("Pre-flight: closed form vs MN on m<=24 (both c) ...")
    pf = 0
    for c in (4, 5):
        for m in range(c + 2, 25):
            n = 2 * m
            for av in range(c + 1, n):
                bv = n - c - av
                if bv <= c or bv >= av:
                    continue
                for jv in range(0, bv + 1):
                    if Mj((av, bv, c), jv, m) != M_closed(c, av, bv, jv, m):
                        pf += 1
    print(f"  pre-flight mismatches = {pf}\n")

    res4 = run_census(4, MAXM)
    res5 = run_census(5, MAXM)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for r in (res4, res5):
        print(f"  c={r['c']}: shapes={r['n_shapes']}  toggles={r['n_tie']}  "
              f"rigid={r['n_shapes']-r['n_tie']}  "
              f"tie-locus={r['n_tie_in_possible']}/{r['n_tie_possible']}  "
              f"MN-vs-closed mism={r['mism']}  box-viol={r['bad_box']}")
