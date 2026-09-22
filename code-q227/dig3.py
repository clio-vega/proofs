from itertools import combinations
from beads import *
n = 3
def all_beadsets(n):
    return [frozenset(c) for k in range(1, n) for c in combinations(range(n), k)]
for S, T, win in additive_pairs(n):
    realA = []
    for A in all_beadsets(n):
        nu = act_word(A, word_cd(T, n), n)
        if nu is None: continue
        lam = act_word(nu, word_cd(S, n), n)
        if lam is not None: realA.append((A, nu, lam))
    if not realA: continue
    nruns = len(runs(S, n))
    mv = clio_moves(S, T, n)
    if len(mv) != nruns:
        print(f"S={sorted(S)} T={sorted(T)} win={win} runs(S)={runs(S,n)} "
              f"letters={sorted(clio_letters(S,T,n))} #mv={len(mv)} nruns={nruns}")
        for (A, nu, lam) in realA:
            print(f"    A={sorted(A)} -> nu={sorted(nu)} -> lam={sorted(lam)}")
