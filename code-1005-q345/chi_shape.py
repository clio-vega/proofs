"""Section 7 of the paper: the shape of chi.  Three claims, each with its check.
(A) closed form at k=2;  (B) shift-invariance of the top two-letter part;
(C) chi is the signed count over (weakly increasing sequence, marked subset)."""
from itertools import combinations_with_replacement, combinations
from math import comb
from chi import h_in_Y
n = 10

badA = totA = 0
for a in range(1, 9):
    P = h_in_Y(a, 2, n)
    for j in range(a+1):
        ex = [0]*(n-1); ex[0] = a-j; ex[1] = j
        pred = sum((-1)**(s-j)*comb(s, j) for s in range(j, a+1))
        totA += 1
        if P.get(tuple(ex), 0) != pred: badA += 1
print(f"(A) k=2 closed form [Y1^(a-j)Y2^j] h_a = sum_{{s=j}}^a (-1)^(s-j) C(s,j): {totA-badA}/{totA}, {badA} bad")

badB = totB = 0
for k in range(2, 6):
    for a in range(1, 7):
        P, Q = h_in_Y(a, k, n), h_in_Y(a, 2, n)
        for j in range(a+1):
            e1 = [0]*(n-1); e1[k-2] = a-j; e1[k-1] = j
            e2 = [0]*(n-1); e2[0] = a-j; e2[1] = j
            totB += 1
            if P.get(tuple(e1), 0) != Q.get(tuple(e2), 0): badB += 1
print(f"(B) shift-invariance of the {{Y_(k-1),Y_k}}-part: {totB-badB}/{totB}, {badB} bad")

badC = totC = signed = 0
for k in range(1, 5):
    for a in range(1, 6):
        acc = {}
        for seq in combinations_with_replacement(range(1, k+1), a):
            for r in range(a+1):
                for S in combinations(range(a), r):
                    ex = [0]*(n-1); ok = True
                    for idx, i in enumerate(seq):
                        v = i - (1 if idx in S else 0)
                        if v == 0: ok = False; break
                        ex[v-1] += 1
                    if ok:
                        signed += 1
                        acc[tuple(ex)] = acc.get(tuple(ex), 0) + (-1)**r
        totC += 1
        if {e: c for e, c in acc.items() if c} != h_in_Y(a, k, n): badC += 1
print(f"(C) signed count over (weakly incr. seq, marked subset): {totC-badC}/{totC}, {badC} bad "
      f"({signed} signed terms collapsing across the 20 cases -- the cancellation is heavy)")

print("\nTable (paper S7):")
for k in range(1, 5):
    for a in range(1, 6):
        P = h_in_Y(a, k, n)
        t = sorted(((tuple(i+1 for i, e in enumerate(ex) for _ in range(e)), c) for ex, c in P.items()),
                   key=lambda z: (len(z[0]), z[0]))
        print(f"  h_{a}(x_1..x_{k}) = " + " ".join(
            f"{'+' if c>0 else '-'}{abs(c) if abs(c)!=1 else ''}Y_{''.join(map(str,w))}" for w, c in t))
