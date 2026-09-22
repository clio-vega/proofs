"""Q227 verdict.  Claims, each stated so that it can FAIL."""
from itertools import combinations
from beads import *
import collections, sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from emptyX import case3_bs, ext_move

def all_beadsets(n):
    return [frozenset(c) for k in range(1, n) for c in combinations(range(n), k)]

def hops(A, n):                       # available forward TASEP hops of a state
    return sum(1 for p in A if (p + 1) % n not in A)

def MS_letters(S, T, n):
    X = X_set(S, T, n); out = set()
    if X:
        for x in X:
            r = ms_etilde(S, T, x, n)
            if r is not None: out.add(r[0])
    else:
        for b, t in case3_bs(S, T, n):
            r = ext_move(S, T, b, n, 'A')
            if r: out.add(r[0])
    return out

def is321(win, n):
    def W(i):
        q, r = divmod(i - 1, n); return win[r] + q * n
    hi = 3 * n
    for x in range(1, n + 1):
        for y in range(x + 1, hi):
            if W(y) >= W(x): continue
            for z in range(y + 1, hi):
                if W(z) < W(y): return False
    return True

hdr = f"{'n':>2} {'pairs':>6} {'real':>6} {'real=321':>9} | {'e=m on locus':>13} {'e>m off':>8} | " \
      f"{'#mv=#usable runs':>17} {'unusable runs have [m+1,M+1]<=T':>32} | " \
      f"{'(C) #mv=#hops(nu)':>18} | {'excess mv':>9} {'excess ON locus':>15}"
print(hdr)
GT = collections.Counter()
for n in range(3, 8):
    beads = all_beadsets(n)
    c = collections.Counter()
    for S, T, win in additive_pairs(n):
        c['pairs'] += 1
        wS, wT = word_cd(S, n), word_cd(T, n)
        realA = []
        for A in beads:
            nu = act_word(A, wT, n)
            if nu is None: continue
            lam = act_word(nu, wS, n)
            if lam is not None: realA.append((A, nu, lam))
        real = bool(realA)
        c['real'] += real
        c['real_eq_321'] += (real == is321(win, n))
        E = clio_letters(S, T, n)
        R = runs(S, n) if len(S) < n else []
        usable = [r for r in R if any(x in E for x in r)]
        e_is_m = all(r[0] in E for r in usable)
        unus_ok = all(all((j + 1) % n in T for j in r) for r in R if r not in usable)
        m_notin_T = all(r[0] not in T for r in R)
        if real:
            c['loc'] += 1
            c['loc_e_is_m'] += e_is_m
            c['loc_m_notin_T'] += m_notin_T
            c['loc_unus_ok'] += unus_ok
            c['loc_mv_eq_usable'] += (len(clio_moves(S, T, n)) == len(usable))
            for (A, nu, lam) in realA:
                c['triples'] += 1
                c['C_hops'] += (len(E) == hops(nu, n))
        else:
            c['off'] += 1
            c['off_e_gt_m'] += (not e_is_m)
        B = MS_letters(S, T, n)
        if E - B:
            c['exc_pairs'] += 1
            c['exc_moves'] += len(E - B)
            if real: c['exc_on_locus'] += 1
    print(f"{n:>2} {c['pairs']:>6} {c['real']:>6} {c['real_eq_321']:>9} | "
          f"{c['loc_e_is_m']}/{c['loc']:<12} {c['off_e_gt_m']:>8} | "
          f"{c['loc_mv_eq_usable']}/{c['loc']:<15} {c['loc_unus_ok']}/{c['loc']:<30} | "
          f"{c['C_hops']}/{c['triples']:<16} | {c['exc_moves']:>9} {c['exc_on_locus']:>15}")
    GT.update(c)
print()
print("TOTALS:", dict(GT))
