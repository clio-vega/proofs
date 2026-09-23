"""Q232: verify the MECHANISM of the proof, not just its conclusion.

The proof asserts four structural facts.  Each is checked here independently
of Shi's length formula, and the final additivity verdict is cross-checked
against Shi's formula (which is what affine.length implements).

  (A)  u_{S'} and u_S differ exactly at positions == p and == i+1 (mod n),
       where [p,q] is the run of S containing i; and
       u_S(p)=q+1, u_{S'}(p)=i, u_S(i+1)=i, u_{S'}(i+1)=q+1.
  (B)  for x <= p-1:  u_S(x) <= p-1  and  u_{S'}(x) <= p-1.
  (P)  p  is NOT congruent to  i+1  mod n.
  (C)  every run [m',M'] of T' falls in exactly one of Case 1a / 1b / 2,
       and in each case the inequality u_{S'}(M'+1) > u_{S'}(j) holds for the
       reason the proof gives.
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import u_S, window, compose, length
from itertools import combinations


def is_additive(S, T, n):
    if len(S) >= n or len(T) >= n:
        return False
    return length(window(compose(u_S(S, n), u_S(T, n)), n), n) == len(S) + len(T)


def zruns(S, n, lo, hi):
    """Maximal integer intervals of the periodic lift of S inside [lo,hi],
    returned only when fully contained (both endpoints have a gap outside)."""
    Sl = set(x % n for x in S)
    out = []
    x = lo
    while x <= hi:
        if x % n in Sl and (x - 1) % n not in Sl:
            m = x
            M = x
            while (M + 1) % n in Sl:
                M += 1
            out.append((m, M))
            x = M + 1
        else:
            x += 1
    return out


def crit(w, T, n):
    """lem:add criterion: for every run [m,M] of T,  w(M+1) > w(j) on [m,M]."""
    if len(T) >= n:
        return None
    for (m, M) in zruns(T, n, 0, n - 1):
        for j in range(m, M + 1):
            if not w(M + 1) > w(j):
                return False
    return True


def subs(n):
    return [frozenset(c) for r in range(n) for c in combinations(range(n), r)]


stats = {'A': 0, 'B': 0, 'P': 0, '1a': 0, '1b': 0, '2': 0,
         'crit_vs_shi': 0, 'fail': 0}
bad = []

for n in range(2, 9):
    sl = subs(n)
    for S in sl:
        for T in sl:
            if not is_additive(S, T, n):
                continue
            uS = u_S(S, n)
            # ---- cross-check lem:add against Shi for the ORIGINAL pair ----
            if crit(uS, T, n) is not True:
                bad.append(('lem:add disagrees with Shi on (S,T)', n, S, T))
            for i in sorted(S):
                Sp = frozenset(S) - {i}
                Tp = frozenset(T) - {(i + 1) % n}
                uSp = u_S(Sp, n)
                # run [p,q] of S containing i, as integers with p <= i <= q
                p = i
                while (p - 1) % n in S:
                    p -= 1
                q = i
                while (q + 1) % n in S:
                    q += 1
                # ---- (P) ----
                if (p - (i + 1)) % n == 0:
                    bad.append(('P fails', n, S, T, i))
                else:
                    stats['P'] += 1
                # ---- (A) ----
                okA = True
                for x in range(p - 2 * n, q + 2 * n + 2):
                    same = (uSp(x) == uS(x))
                    differs_expected = (x % n == p % n) or (x % n == (i + 1) % n)
                    if same == differs_expected:
                        okA = False
                if not (uS(p) == q + 1 and uSp(p) == i
                        and uS(i + 1) == i and uSp(i + 1) == q + 1):
                    okA = False
                if okA:
                    stats['A'] += 1
                else:
                    bad.append(('A fails', n, S, T, i))
                # ---- (B) ----
                okB = all(uS(x) <= p - 1 and uSp(x) <= p - 1
                          for x in range(p - 2 * n, p))
                if okB:
                    stats['B'] += 1
                else:
                    bad.append(('B fails', n, S, T, i))
                # ---- (C): walk every run of T' and classify ----
                if len(Tp) < n:
                    for (m2, M2) in zruns(Tp, n, -2 * n, 3 * n):
                        if not (0 <= m2 <= n - 1):
                            continue                     # one lift per residue
                        if (M2 + 1) % n == p % n:
                            case = '2'
                        elif (M2 + 1) % n not in set(x % n for x in T):
                            case = '1a'
                        else:
                            case = '1b'
                            if (M2 + 1) % n != (i + 1) % n:
                                bad.append(('case split incomplete', n, S, T, i))
                        stats[case] += 1
                        for j in range(m2, M2 + 1):
                            if not uSp(M2 + 1) > uSp(j):
                                bad.append(('criterion fails', n, S, T, i, case,
                                            m2, M2, j))
                                stats['fail'] += 1
                        # case-specific reason the proof gives
                        if case == '2':
                            sh = (M2 + 1) - p            # multiple of n
                            if not all(uSp(j) <= p - 1 + sh for j in range(m2, M2 + 1)):
                                bad.append(('case2 reason fails', n, S, T, i))
                        if case == '1a':
                            if not all(uSp(j) <= uS(j) for j in range(m2, M2 + 1)):
                                bad.append(('1a monotone-down fails', n, S, T, i))
                            if not uSp(M2 + 1) >= uS(M2 + 1):
                                bad.append(('1a monotone-up fails', n, S, T, i))
                # ---- lem:add vs Shi for the REDUCED pair ----
                c = crit(uSp, Tp, n)
                s = is_additive(Sp, Tp, n)
                if c != s:
                    bad.append(('lem:add disagrees with Shi on (S-i,T-(i+1))',
                                n, S, T, i))
                else:
                    stats['crit_vs_shi'] += 1
                if not s:
                    bad.append(('THEOREM FALSE', n, sorted(S), sorted(T), i))

print("structural checks (all additive (S,T) and all i in S, n=2..8):")
for k in ('P', 'A', 'B', '1a', '1b', '2', 'crit_vs_shi', 'fail'):
    print(f"  {k:>12}: {stats[k]}")
print(f"violations: {len(bad)}")
for b in bad[:10]:
    print("  ", b)
