"""Q254, Tests A and B: does the NATURAL column dictionary work?

Test A (calibration, d=0).  For an ordinary skew shape lam/mu, the columns
D_j = {i : mu_i < j <= lam_i} are a diagram in FMS's sense.  Is
supp(chi_D) = supp(s_{lam/mu}(x_1..x_ell))?

Test B (the real one, d>=1).  Same question for a cylindric skew shape, whose
columns are read off one fundamental domain.

Negative controls: drop one column; reverse the column order.
"""
import sys, itertools
from collections import defaultdict
sys.path.insert(0, '/home/clio/projects/proofs/code-q254')
from fms import supp_chi, active_rows, is_Sm_stable, is_bottom_justified
import cyl as C   # vendored copy, see cyl.py header


# ---------- ordinary skew Schur support, computed independently -------------
def ssyt_weights(lam, mu, ell):
    """supp(s_{lam/mu}(x_1..x_ell)) by direct SSYT enumeration."""
    lam = list(lam); mu = list(mu) + [0]*(len(lam)-len(mu))
    rows = len(lam)
    cells = [(i, j) for i in range(rows) for j in range(mu[i], lam[i])]
    out = set()
    def rec(k, T):
        if k == len(cells):
            w = [0]*ell
            for v in T.values(): w[v-1] += 1
            out.add(tuple(w)); return
        i, j = cells[k]
        lo = 1
        if j > mu[i] and (i, j-1) in T: lo = max(lo, T[(i, j-1)])      # weakly incr in row
        if i > 0 and (i-1, j) in T: lo = max(lo, T[(i-1, j)] + 1)      # strictly incr in col
        for v in range(lo, ell+1):
            T[(i, j)] = v; rec(k+1, T); del T[(i, j)]
    rec(0, {})
    return out


def skew_columns(lam, mu, shift=0):
    lam = list(lam); mu = list(mu) + [0]*(len(lam)-len(mu))
    return [tuple(i+1+shift for i in range(len(lam)) if mu[i] < j <= lam[i])
            for j in range(1, (lam[0] if lam else 0) + 1)]


def partitions_upto(N, maxparts):
    out = []
    def rec(rem, mx, cur):
        out.append(tuple(cur))
        for p in range(min(rem, mx), 0, -1):
            if len(cur) < maxparts: rec(rem-p, p, cur+[p])
    rec(N, N, [])
    return sorted(set(out))


def test_A(maxsize=8, maxrows=4):
    tot = agree = sym = 0
    first = None
    for rows in range(1, maxrows+1):
        for lam in partitions_upto(maxsize, rows):
            if len(lam) != rows: continue
            for mu in partitions_upto(sum(lam), rows):
                if len(mu) > rows: continue
                mup = list(mu)+[0]*(rows-len(mu))
                if any(mup[i] > lam[i] for i in range(rows)): continue
                if not all(mup[i] >= mup[i+1] for i in range(rows-1)): continue
                N = sum(lam)-sum(mup)
                if N == 0 or N > maxsize: continue
                ell = rows
                D = skew_columns(lam, mup)
                if active_rows(D) != ell: continue
                tot += 1
                S_fms = supp_chi(D, ell)
                S_true = ssyt_weights(lam, mup, ell)
                if S_fms == S_true: agree += 1
                elif first is None: first = (lam, tuple(mup), ell, D, sorted(S_fms), sorted(S_true))
                if is_Sm_stable(S_fms, ell): sym += 1
    print(f"TEST A (d=0, natural column diagram): {agree}/{tot} shapes where "
          f"supp(chi_D) = supp(s_lam/mu);  {sym}/{tot} have S_ell-stable supp(chi_D)")
    if first:
        lam, mu, ell, D, a, b = first
        print(f"   minimal witness: lam={lam}, mu={mu}, ell={ell}, columns={D}")
        print(f"      supp(chi_D)      = {a}")
        print(f"      supp(s_lam/mu)   = {b}")
    return tot, agree, sym


def test_A_all_shifts(maxsize=6, maxrows=3, maxell=5):
    """Stronger: allow the skew shape to be placed ANYWHERE in the grid."""
    tot = hit = 0
    for rows in range(1, maxrows+1):
        for lam in partitions_upto(maxsize, rows):
            if len(lam) != rows: continue
            for mu in partitions_upto(sum(lam), rows):
                mup = list(mu)+[0]*(rows-len(mu))
                if len(mu) > rows or any(mup[i] > lam[i] for i in range(rows)): continue
                if not all(mup[i] >= mup[i+1] for i in range(rows-1)): continue
                N = sum(lam)-sum(mup)
                if N == 0 or N > maxsize: continue
                for ell in range(rows, maxell+1):
                    S_true = ssyt_weights(lam, mup, ell)
                    if not S_true: continue
                    tot += 1
                    ok = False
                    for shift in range(0, ell-rows+1):
                        D = skew_columns(lam, mup, shift)
                        if supp_chi(D, ell) == S_true: ok = True; break
                    if ok: hit += 1
    print(f"TEST A' (d=0, ANY vertical placement of the skew diagram): "
          f"{hit}/{tot} (lam/mu,ell) matched by some shift")
    return tot, hit


if __name__ == "__main__":
    test_A()
    test_A_all_shifts()
