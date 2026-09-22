"""Q220 STEP 1: compare  U_{x in X} {e~_1^{(x)}(uv)}  with  {Thm 4.1 moves on uv}.

Both sides are computed as sets of pairs (S', T').  We also record the
"letter sets": the letter removed from S on each side.
"""
import sys, collections
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import *

def report(n, verbose_limit=6):
    stats = collections.Counter()
    excess_clio = []    # Thm 4.1 moves that are no e~_1^{(x)} for any x
    excess_ms   = []    # e~_1^{(x)} moves that are not Thm 4.1 moves
    noX         = []    # X empty but Thm 4.1 has a move
    for S, T, wv in additive_pairs(n):
        X = X_set(S, T, n)
        E = clio_moves(S, T, n)
        B = set()
        for x in X:
            r = ms_etilde(S, T, x, n)
            if r is not None:
                B.add((r[1], r[2]))
        stats['pairs'] += 1
        if not X:
            stats['X_empty'] += 1
            if E:
                noX.append((S, T, E))
        for mv in E - B:
            excess_clio.append((S, T, mv, X))
        for mv in B - E:
            excess_ms.append((S, T, mv, X))
        if E == B:
            stats['equal'] += 1
        stats['E_moves'] += len(E)
        stats['B_moves'] += len(B)
    return stats, excess_clio, excess_ms, noX

for n in range(3, 8):
    stats, ec, em, noX = report(n)
    print(f"n={n}: additive pairs {stats['pairs']}, X empty {stats['X_empty']}, "
          f"E==B on {stats['equal']}, |E| total {stats['E_moves']}, |B| total {stats['B_moves']}")
    print(f"      Thm4.1-not-MS: {len(ec)}   MS-not-Thm4.1: {len(em)}   (X empty but E nonempty: {len(noX)})")
    for S,T,mv,X in ec[:4]:
        print(f"        [E\\B] S={sorted(S)} T={sorted(T)} X={X} -> S'={sorted(mv[0])} T'={sorted(mv[1])}")
    for S,T,mv,X in em[:4]:
        print(f"        [B\\E] S={sorted(S)} T={sorted(T)} X={X} -> S'={sorted(mv[0])} T'={sorted(mv[1])}")
