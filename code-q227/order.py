"""Which factor of v = u_S u_T acts on the beads FIRST?

Convention check, done against the nilCoxeter identity, not fitted:
u_S u_T = u_{S'} u_{T'} = u_v when both factorisations are length-additive
(Theorem 4.1).  So the two words must have the SAME action on every A.
"""
from itertools import combinations
from beads import *

def all_beadsets(n):
    out = []
    for k in range(1, n):
        for c in combinations(range(n), k):
            out.append(frozenset(c))
    return out

for n in range(3, 7):
    tot = ok = real = 0
    firstS = firstT = 0
    for S, T, w in additive_pairs(n):
        mvs = clio_moves(S, T, n)
        if not mvs:
            continue
        word_v = word_cd(S, n) + word_cd(T, n)
        for A in all_beadsets(n):
            img = act_word(A, word_v, n)
            if img is None:
                continue
            real += 1
            # does u_T act first (rightmost-first on the concatenated word)?
            midT = act_cd(A, T, n)
            midS = act_cd(A, S, n)
            if midT is not None and act_cd(midT, S, n) == img:
                firstT += 1
            if midS is not None and act_cd(midS, T, n) == img:
                firstS += 1
            for (Sp, Tp) in mvs:
                tot += 1
                img2 = act_word(A, word_cd(Sp, n) + word_cd(Tp, n), n)
                ok += (img2 == img)
    print(f"n={n}: realisable (S,T,A) triples with >=1 move: {real};"
          f"  T-acts-first {firstT}, S-acts-first {firstS};"
          f"  exchange-move image agrees: {ok}/{tot}")
