"""PROVE 2026-11-13 c2, Step 0 -- the pre-registered w0/evacuation falsifier for the DREAM connection
'the backward pairing is the longest element'.

Objects (conventions of code-1113-Q585-H1/trace2.py):
  T      : beta_CM LR tableau of skew shape nu/mu, content lambda; cells (c,r), r = row index.
  journeys g[0..m-1] : the G-phase migration journeys (mu-nest), in the engine's order.
  slides  : Sigma = switching of (U_mu , T), U_mu the standard tableau of shape mu giving the slide order:
            the slide into the cell carrying the LARGEST entry of the order-tableau is performed first.
  U      : row-superstandard SYT of shape mu (1..mu_1 along row 0, ...).  This is the order sigma_moves used.

Pairings tested:
  CTRL-i  (known failure): order U, journey k <-> slide m+1-k.  Must reproduce 'R': 20 step-strings.
  P       (pre-registered): order xi(U), journey k <-> slide k.
  CTRL-ii (wrong factor): order U, T replaced by its Schutzenberger (Lusztig) involution xi_cryst(T) (same skew
           shape, lowest weight of T's crystal component), journey k <-> slide k.  Must differ from P somewhere.
  CTRL-iii: xi(xi(U)) == U and xi(U) is an SYT of shape mu; xi_cryst o xi_cryst = id, content reversed.
Falsified iff under P some move-string carries >= 2 step-strings.
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0, '/home/clio/projects/proofs/code-1113-Q585-H1')
from trace2 import gtrace, parts, LR, beta_cm, trim
from collections import defaultdict, Counter

# ---------------------------------------------------------------- straight-shape evacuation
def shape_rows(cells):
    rows = Counter(r for c, r in cells); return [rows[r] for r in range(len(rows))]

def row_superstandard(mu):
    U = {}; k = 0
    for r, l in enumerate(mu):
        for c in range(l): k += 1; U[(c, r)] = k
    return U

def evacuate(S):
    """Schutzenberger evacuation of an SYT S (dict (c,r)->entry, straight shape, English, r down).
    Repeatedly delete the entry 1 (at (0,0)), jdt-slide the hole outward (smaller neighbour moves in),
    subtract 1 from all entries; the cell finally vacated at step t receives n+1-t."""
    S = dict(S); n = len(S); out = {}
    for t in range(1, n + 1):
        c, r = 0, 0; del S[(0, 0)]
        while True:
            R = S.get((c + 1, r)); D = S.get((c, r + 1))
            if R is None and D is None: break
            if D is not None and (R is None or D < R):
                S[(c, r)] = D; del S[(c, r + 1)]; r += 1
            else:
                S[(c, r)] = R; del S[(c + 1, r)]; c += 1
        out[(c, r)] = n + 1 - t
        S = {x: v - 1 for x, v in S.items()}
    return out

def is_syt(S):
    cells = set(S)
    for (c, r), v in S.items():
        if (c + 1, r) in cells and S[(c + 1, r)] <= v: return False
        if (c, r + 1) in cells and S[(c, r + 1)] <= v: return False
        if c > 0 and (c - 1, r) not in cells: return False
        if r > 0 and (c, r - 1) not in cells: return False
    return sorted(S.values()) == list(range(1, len(S) + 1))

# ---------------------------------------------------------------- switching with a given slide order
def sigma_moves_order(S, outer, order_tab):
    """Same jdt as trace2.sigma_moves, but slides go into the cells of mu in DEcreasing order of order_tab."""
    S = dict(S); res = []
    cells = sorted(order_tab, key=lambda x: -order_tab[x])
    for (c, r) in cells:
        hc, hr = c, r; mv = ''
        while True:
            R = S.get((hc + 1, hr)); D = S.get((hc, hr + 1))
            if R is None and D is None: break
            if D is not None and (R is None or D <= R):
                S[(hc, hr)] = D; del S[(hc, hr + 1)]; hr += 1; mv += 'U'
            else:
                S[(hc, hr)] = R; del S[(hc + 1, hr)]; hc += 1; mv += 'R'
        res.append(((c, r), (hc, hr), mv))
    return res, S

# ---------------------------------------------------------------- crystal Schutzenberger involution on skew T
def word_cells(T):   # same reading word as test_K.py: rows in decreasing r, each left->right
    rows = sorted({r for c, r in T}, reverse=True)
    return [(c, r) for r in rows for c in sorted(c for c, rr in T if rr == r)]

def f_i(T, i):
    """lower: rightmost unmatched i becomes i+1 (bracketing (i+1 ... i) as in test_K.e_i)."""
    cells = word_cells(T); unmatched_i = []
    # scan right-to-left: an i is matched by an i+1 to its LEFT ... use standard rule mirrored from e_i:
    # in test_K, scanning left->right, i+1 pushed, i pops: pairs are (i+1 left, i right). Unmatched i's are
    # those with no unmatched i+1 to their left; f_i changes the RIGHTMOST?  No: the LEFTMOST? -- by the
    # signature rule the unmatched subword is i^a (i+1)^b?  Here the pairing (i+1, i) leaves unmatched word
    # i^a (i+1)^b; e_i changes the leftmost unmatched i+1 and f_i the rightmost unmatched i.
    stack = []; unm_i = []
    for x in cells:
        v = T[x]
        if v == i + 1: stack.append(x)
        elif v == i:
            if stack: stack.pop()
            else: unm_i.append(x)
    if not unm_i: return None
    T = dict(T); T[unm_i[-1]] = i + 1; return T

def e_i(T, i):
    cells = word_cells(T); stack = []
    for x in cells:
        v = T[x]
        if v == i + 1: stack.append(x)
        elif v == i and stack: stack.pop()
    if not stack: return None
    T = dict(T); T[stack[0]] = i; return T

def lowest(T, n):
    changed = True
    while changed:
        changed = False
        for i in range(1, n):
            while True:
                T2 = f_i(T, i)
                if T2 is None: break
                T = T2; changed = True
    return T

def highest(T, n):
    changed = True
    while changed:
        changed = False
        for i in range(1, n):
            while True:
                T2 = e_i(T, i)
                if T2 is None: break
                T = T2; changed = True
    return T

def is_ssyt(T):
    for (c, r), v in T.items():
        if (c + 1, r) in T and T[(c + 1, r)] < v: return False
        if (c, r + 1) in T and T[(c, r + 1)] <= v: return False
    return True

# ---------------------------------------------------------------- census
def census(amax, bmax, cmin, gl=4):
    for a in range(1, amax + 1):
        for b_ in range(1, bmax + 1):
            for lam in parts(a, gl):
                for mu in parts(b_, gl):
                    for nu, bs in LR(lam, mu, gl).items():
                        if len(bs) < cmin: continue
                        nt = trim(nu)
                        for bb in bs: yield lam, mu, nt, beta_cm(bb, lam, mu, nt)

def run(amax, bmax, cmin, label):
    tab = {k: defaultdict(set) for k in ('ctrl_i', 'P', 'ctrl_ii', 'P_rev')}
    land = Counter(); n = 0; diff_ii = 0; xi_nontriv = 0; chk = Counter(); hits = {k: Counter() for k in tab}
    for lam, mu, nt, T in census(amax, bmax, cmin):
        g = gtrace(T, nt); m = len(g)
        mu_t = tuple(x for x in mu if x)
        U = row_superstandard(mu_t); X = evacuate(U)
        chk['xi syt'] += is_syt(X); chk['xi^2=id'] += (evacuate(X) == U); chk['U'] += 1
        xi_nontriv += (X != U)
        sU, _ = sigma_moves_order(T, nt, U); sX, _ = sigma_moves_order(T, nt, X)
        l = len(lam); Tl = lowest(T, max(l, max(T.values())))
        chk['xc ssyt'] += is_ssyt(Tl); chk['xc^2'] += (highest(Tl, l) == T)
        chk['xc content reversed'] += (sorted(Counter(Tl.values()).items()) ==
                                       sorted((l + 1 - v, k) for v, k in Counter(T.values()).items()))
        sL, _ = sigma_moves_order(Tl, nt, U)
        assert len(sU) == m == len(sX) == len(sL), (lam, mu, nt, m, len(sU))
        for k in range(m):
            v, s, cell = g[k]; n += 1
            for key, mv in (('ctrl_i', sU[m - 1 - k][2]), ('P', sX[k][2]), ('ctrl_ii', sL[k][2]), ('P_rev', sX[m - 1 - k][2])):
                tab[key][mv].add(s); hits[key][mv] += 1
            land['ctrl_i'] += (cell == sU[m - 1 - k][1]); land['P'] += (cell == sX[k][1])
            land['P_rev'] += (cell == sX[m - 1 - k][1])
            diff_ii += (sX[k][2] != sL[k][2])
    print('==', label, 'journeys', n, 'xi(U)!=U on', xi_nontriv, 'elements; checks', dict(chk))
    for key, d in tab.items():
        bad = {mv: len(ss) for mv, ss in d.items() if len(ss) > 1}
        rep = sum(1 for mv in d if hits[key][mv] >= 2)
        print(f'  {key:8s}: distinct move-strings {len(d)}, reached>=2x {rep}, move-strings with >=2 step-strings {len(bad)}',
              ' detail', dict(sorted({mv: len(ss) for mv, ss in d.items()}.items())))
    print('  landing cell matches:', dict(land), ' ctrl_ii differs from P on', diff_ii, 'journeys')
    return tab

if __name__ == '__main__':
    run(4, 4, 2, 'c>=2, |lam|,|mu|<=4 (the 74)')
    run(5, 5, 3, 'c>=3, |lam|,|mu|<=5')
