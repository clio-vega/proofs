"""Direct test of Lemma 1: u_S . A != 0  iff  for every run [m,M] of S,
   m in A and [m+1,M+1] cap A = empty; and then the bead at m travels to M+1."""
from itertools import combinations
from beads import *
ok = bad = 0
for n in range(3, 9):
    subs = [frozenset(c) for k in range(0, n) for c in combinations(range(n), k)]
    for S in subs:
        for A in [frozenset(c) for k in range(1, n) for c in combinations(range(n), k)]:
            pred_ok = all(r[0] in A and all((j) % n not in A
                          for j in range(r[0] + 1, r[0] + len(r) + 1)) for r in runs(S, n))
            img = act_cd(A, S, n)
            if (img is not None) != pred_ok:
                bad += 1; continue
            if img is not None:
                want = set(A)
                for r in runs(S, n):
                    want.discard(r[0]); want.add((r[0] + len(r)) % n)
                if frozenset(want) != img:
                    bad += 1; continue
            ok += 1
print(f"Lemma 1 checked on {ok+bad} (S,A) pairs, n=3..8:  agree {ok}, disagree {bad}")
