"""(G1'): transporting the reduced-pair exchange move back to (S,T).

Setting: (S,T) additive, S u T = Z/n (so X(S,T) = {}), b a Case-(3) index:
    [b-t, b-1] subset S,  b not in S,  [b-t+1, b] subset T,  b-t not in T,  t>=1.
MS reduce to  So = S-{b-1},  To = T-{b}  and proceed with cut x=b.
Let c = min_{<_b} L_1 computed on (So,To), and m = bottom of the run of c.

CLAIMED PROOF (checked step by step below):
 (0) readings A and B always agree, because b-1 is the TOP of its S-run, so
     deleting it changes no run's BOTTOM.
 (1) b is in X(So,To), and (So,To) is additive by the deletion theorem, so
     Thm 3.2 applies to (So,To) at x=b and gives
        c+1 not in To,   [mo+1,c] subset To,   mo not in To.
 (2) transport: m = mo;  m not in T (since m not in To and m != b, as m in S);
     c+1 not in T (since c+1 not in To and c+1 != b, as c != b-1);
     every j in [m,c-1] has j+1 in To subset T.
     The run of c in S is the run of c in So, except possibly extended at the
     TOP from b-2 to b-1 -- which cannot affect m or the FIRST index.
     Hence c = e_{[m,M]} for the run [m,M] of c in S, that run is usable,
     and (S-{c}, T+{m}) is the Thm 4.1 move on it.
 (3) side condition of Lem 4.3 for (S,T):  |T| <= n-2.
"""
import sys, collections
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import (u_S, window, compose, length, runs, clio_letters,
                    clio_moves, ms_order, ms_pairing, X_set, additive_pairs)
from emptyX import case3_bs


def run_of(e, S, n):
    for r in runs(S, n):
        if e % n in r:
            return r
    return None


c = collections.Counter()
bad = collections.defaultdict(list)


def note(tag, *info):
    c[tag] += 1
    if len(bad[tag]) < 4:
        bad[tag].append(info)


for n in range(3, 9):
    for S, T, _ in additive_pairs(n):
        if X_set(S, T, n):
            continue
        c['pairs'] += 1
        for b, t in case3_bs(S, T, n):
            c['case3'] += 1
            So = frozenset(S) - {(b - 1) % n}
            To = frozenset(T) - {b % n}
            # --- (1a) reduced pair additive (the theorem proved today) ---
            if length(window(compose(u_S(So, n), u_S(To, n)), n), n) != len(So) + len(To):
                note('reduced-not-additive', n, sorted(S), sorted(T), b)
                continue
            c['reduced-additive'] += 1
            # --- (1b) b in X(So,To) ---
            if b % n in So or b % n in To:
                note('b-not-admissible', n, sorted(S), sorted(T), b)
            else:
                c['b-admissible'] += 1
            # --- (3) side condition |T| <= n-2 ---
            if len(T) > n - 2:
                note('T-too-big', n, sorted(S), sorted(T), b, len(T))
            else:
                c['T-small'] += 1
            L1, _, _ = ms_pairing(So, To, b % n, n)
            if not L1:
                c['L1-empty'] += 1
                continue
            c['L1-nonempty'] += 1
            pos = {r: i for i, r in enumerate(ms_order(b % n, n))}
            cc = min(L1, key=lambda r: pos[r])
            # --- (0) readings agree ---
            rS = run_of(cc, S, n)
            rSo = run_of(cc, So, n)
            if rS is None or rSo is None:
                note('no-run', n, sorted(S), sorted(T), b, cc)
                continue
            mA, mB = rS[0], rSo[0]
            if mA != mB:
                note('readings-differ', n, sorted(S), sorted(T), b, cc, mA, mB)
            else:
                c['readings-agree'] += 1
            # --- (1) Thm 3.2 conclusions on the REDUCED pair ---
            mo, Mo = rSo[0], rSo[-1]
            ok = ((cc + 1) % n not in To and mo not in To
                  and all((j % n) in To for j in range(
                      mo + 1, mo + 1 + ((cc - mo) % n))))
            if ok:
                c['thm32-holds'] += 1
            else:
                note('thm32-fails', n, sorted(S), sorted(T), b, cc)
            # --- (2) the transported claims on the ORIGINAL pair ---
            if mA in set(T):
                note('m-in-T', n, sorted(S), sorted(T), b, cc, mA)
            if (cc + 1) % n in set(T):
                note('c+1-in-T', n, sorted(S), sorted(T), b, cc)
            # run of c in S extends run in So only at the top?
            if not (mA == mB and (rS[-1] == rSo[-1] or rS[-1] == (b - 1) % n)):
                note('run-extension-not-at-top', n, sorted(S), sorted(T), b, cc,
                     rS, rSo)
            else:
                c['run-extends-at-top-only'] += 1
            # --- the conclusion: c is a Thm 4.1 letter and the move is a Thm 4.1 move ---
            if cc not in clio_letters(S, T, n):
                note('c-not-in-E', n, sorted(S), sorted(T), b, cc)
            else:
                c['c-in-E'] += 1
            mv = (frozenset(S) - {cc}, frozenset(T) | {mA})
            if mv not in clio_moves(S, T, n):
                note('move-not-in-Emv', n, sorted(S), sorted(T), b, cc, mA)
            else:
                c['move-in-Emv'] += 1

print("(G1') step-by-step check, n=3..8, all additive (S,T) with X empty:")
for k in ('pairs', 'case3', 'reduced-additive', 'b-admissible', 'T-small',
          'L1-empty', 'L1-nonempty', 'readings-agree', 'thm32-holds',
          'run-extends-at-top-only', 'c-in-E', 'move-in-Emv'):
    print(f"  {k:>26}: {c[k]}")
print("  failures:")
tot = 0
for k, v in bad.items():
    print(f"    {k}: {c[k]}")
    for x in v:
        print("       ", x)
    tot += c[k]
print(f"  total failures: {tot}")


# ---------------------------------------------------------------------------
# The one real hole found above:  |T| <= n-2 FAILS in 894 of 6786 Case-(3)
# instances.  Claim:  in every such instance  L_1 is EMPTY, so no move is
# produced and nothing needs transporting.  Proof (checked here):
#   |T| = n-1 gives T_o = T - {b} exactly two gaps, b and g (the gap of T).
#   Thm 3.2 would give  m_o not in T_o  and  c+1 not in T_o, so each of
#   m_o, c+1 lies in {b, g}.  m_o != b since m_o in S and b not in S, so
#   m_o = g.  And c+1 != b since c != b-1 (c lies in S_o).  So c+1 = g = m_o,
#   i.e. c = m_o - 1 -- impossible, since c lies in the run [m_o, M_o] whose
#   length is at most n-1.
print()
print("the |T| = n-1 instances:")
d = collections.Counter()
for n in range(3, 9):
    for S, T, _ in additive_pairs(n):
        if X_set(S, T, n):
            continue
        for b, t in case3_bs(S, T, n):
            if len(T) <= n - 2:
                continue
            d['T=n-1'] += 1
            So = frozenset(S) - {(b - 1) % n}
            To = frozenset(T) - {b % n}
            L1, _, _ = ms_pairing(So, To, b % n, n)
            d['L1 empty' if not L1 else 'L1 NONEMPTY -- HOLE'] += 1
            gaps = [z for z in range(n) if z not in To]
            if sorted(gaps) != sorted({b % n, [z for z in range(n)
                                               if z not in set(T)][0]}):
                d['gap-structure-wrong'] += 1
            if clio_letters(S, T, n):
                d['E(S,T) NONEMPTY -- HOLE'] += 1
            else:
                d['E(S,T) empty'] += 1
for k in sorted(d):
    print(f"  {k:>26}: {d[k]}")
