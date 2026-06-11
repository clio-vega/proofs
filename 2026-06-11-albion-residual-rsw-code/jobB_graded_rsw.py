"""
JOB B (2026-06-11) -- the GRADED object  G_lambda(q) = sum_{T in SYT(lambda)} q^{s(T)}
                      vs the principal specialisation  s_lambda(1,q,...,q^{k-1}).

This differs from the already-run jobC_rsw (which built G via the symmetric-function
master identity <s_lam, psi^m>).  Here G_lambda(q) is built DIRECTLY from the
parity-twisted SYT descent statistic
    s(T) = sum_{i in Des(T)} w_i,    w_i = 2i-1 if (n-i) odd else 0,
the "branch-exponent" statistic ({d_j} = {s(T)}).

Steps:
  0. CROSS-CHECK (new): does the graded sum_T q^{s(T)} equal the master-identity
     polynomial G_poly(lambda)(t)?  (verifies the two representations are ONE object).
  1. Principal specialisation s_lambda(1,q,...,q^{N-1}) for a menu of N.
  2. At q = zeta_d, d = 2,3,4,5,6: is  G_lambda(zeta_d) = zeta_d^{c} * s_lambda(1,zeta_d,...) ?
     (uniform q-shift c).  Report exact match / match-up-to-shift / discrepancy.
  3. The d=4 vanisher (2,2) specifically.
"""
import sys, csv, os
import sympy as sp
from functools import lru_cache
from collections import defaultdict
from fiber_engine import G_poly, princ_spec_poly, partitions, t as TSYM

RESULTS = "/home/clio/projects/code/results"
q = sp.Symbol('q')

# ----------------------------------------------------------------- SYT machinery
def syt_list(lam):
    """All standard Young tableaux of shape lam, as dict cell->value via filling."""
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    cells = [(i, j) for i, r in enumerate(lam) for j in range(r)]
    results = []
    # build by placing 1..n; a cell is fillable if left and up neighbors filled
    filled = {}
    def rows_cols():
        return lam
    def backtrack(v):
        if v > n:
            results.append(dict(filled))
            return
        for (i, j) in cells:
            if (i, j) in filled:
                continue
            if j > 0 and (i, j-1) not in filled:
                continue
            if i > 0 and (i-1, j) not in filled:
                continue
            filled[(i, j)] = v
            backtrack(v + 1)
            del filled[(i, j)]
    backtrack(1)
    return results

def descents(T):
    """Des(T) = { i : i+1 sits in a strictly lower row than i }."""
    pos = {v: (i, j) for (i, j), v in T.items()}
    n = len(T)
    return [i for i in range(1, n) if pos[i+1][0] > pos[i][0]]

def s_stat(T, n):
    s = 0
    for i in descents(T):
        if (n - i) % 2 == 1:
            s += 2*i - 1
    return s

@lru_cache(maxsize=None)
def graded_G(lam):
    """sum_{T} q^{s(T)}  as an exact sympy poly in q."""
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    poly = sp.Integer(0)
    for T in syt_list(lam):
        poly += q**s_stat(T, n)
    return sp.expand(poly)

# ----------------------------------------------------------------- roots of unity
def zeta_exact(d):
    return sp.exp(2*sp.pi*sp.I/d)

def cval(expr_in_var, var, d):
    """exact value of poly(expr) at var = primitive d-th root of unity, simplified."""
    return sp.simplify(expr_in_var.subs(var, zeta_exact(d)))

def main(NMAX):
    # ---------- STEP 0: representation cross-check ----------
    print("="*72)
    print("STEP 0  graded sum_T q^{s(T)}  ==  master-identity G_poly(lambda)(t) ?")
    print("="*72)
    ok = bad = 0
    mism = []
    for n in range(1, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            gg = graded_G(lam)                       # in q
            mm = sp.expand(G_poly(lam))              # in t
            diff = sp.expand(gg.subs(q, TSYM) - mm)
            if diff == 0:
                ok += 1
            else:
                bad += 1
                if len(mism) < 12:
                    mism.append((lam, sp.expand(gg.subs(q, TSYM)), mm))
    print(f"  match: {ok}   mismatch: {bad}   (lambda |- n, n<={NMAX})")
    for lam, g, m in mism:
        print(f"    MISMATCH lam={lam}:  graded={g}   master={m}")
    graded_eq_master = (bad == 0)

    # ---------- STEP 1+2: principal specialisation comparison ----------
    print("\n" + "="*72)
    print("STEP 1-2  G_lambda(zeta_d) =? zeta_d^c * s_lambda(1,zeta_d,...,zeta_d^{N-1})")
    print("="*72)
    rows = []
    # candidate variable counts N for the principal specialisation
    def Nmenu(lam):
        lam = tuple(x for x in lam if x>0)
        return sorted(set([len(lam), len(lam)+1, lam[0], sum(lam)]))
    ds = [2, 3, 4, 5, 6]
    hit = defaultdict(int); seen = defaultdict(int); hit0 = defaultdict(int)
    for n in range(1, NMAX+1):
        for lam in partitions(n):
            lam = tuple(lam)
            G = graded_G(lam)
            for d in ds:
                zc = zeta_exact(d)
                Lval = sp.simplify(G.subs(q, zc))
                Lzero = (Lval == 0)
                for N in Nmenu(lam):
                    ps = princ_spec_poly(lam, N)
                    if ps == 0:
                        Rval = sp.Integer(0)
                    else:
                        poly, qv = ps
                        Rval = sp.simplify(poly.subs(qv, zc))
                    Rzero = (Rval == 0)
                    seen[(d, N)] += 1
                    best_c = None
                    if Lzero and Rzero:
                        best_c = 0
                    elif (not Lzero) and (not Rzero):
                        for c in range(d):
                            if sp.simplify(Lval - (zc**c)*Rval) == 0:
                                best_c = c; break
                    if best_c is not None:
                        hit[(d, N)] += 1
                        if best_c == 0:
                            hit0[(d, N)] += 1
                    rows.append(dict(n=n, lam='|'.join(map(str,lam)), d=d, N=N,
                                     Lzero=int(Lzero), Rzero=int(Rzero),
                                     shift=('' if best_c is None else best_c)))
    # report: best N per d
    print(f"  {'d':>2} {'bestN':>6} {'hit/seen(any shift)':>20} {'hit/seen(shift0)':>18}")
    for d in ds:
        Ns = sorted(set(N for (dd,N) in seen if dd==d))
        best = max(Ns, key=lambda N: hit[(d,N)])
        print(f"  {d:>2} {best:>6} {f'{hit[(d,best)]}/{seen[(d,best)]}':>20} "
              f"{f'{hit0[(d,best)]}/{seen[(d,best)]}':>18}")

    # ---------- STEP 3: the (2,2) d=4 vanisher ----------
    print("\n" + "="*72)
    print("STEP 3  the d=4 vanisher (2,2)")
    print("="*72)
    lam = (2,2)
    G = graded_G(lam)
    print(f"  graded G_(2,2)(q) = {G}")
    for d in ds:
        print(f"    G_(2,2)(zeta_{d}) = {sp.simplify(G.subs(q, zeta_exact(d)))}")
    for N in Nmenu(lam):
        ps = princ_spec_poly(lam, N)
        if ps != 0:
            poly, qv = ps
            print(f"    s_(2,2)(1,q,..q^{N-1}) = {sp.expand(poly)}")
            print(f"        at zeta_4 = {sp.simplify(poly.subs(qv, zeta_exact(4)))}")

    os.makedirs(RESULTS, exist_ok=True)
    with open(f"{RESULTS}/jobB-graded-rsw.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  wrote results/jobB-graded-rsw.csv ({len(rows)} rows)")
    print(f"\n  STEP 0 verdict: graded == master identity : {graded_eq_master}")

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    main(NMAX)
