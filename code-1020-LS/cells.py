r"""Gate 0b of the PROVE 2026-10-20 brief: SPLIT the hypothesis class of hyp:L^S
BEFORE measuring it.

DOMAIN (stated before any number):
  X in {B,C};  N = 1..5;  every m in [N];  every admissible eps (eps^h = 1,
  h = 2 for C so eps in {+1,-1}; h = 1 for B so eps = +1);
  every ordered pair (A,B) of subsets of [N] with A >= B and m not in A.
  Combinatorics only -- no weight function is evaluated here.

CELLS (the four sub-families the brief says I have never separated):
  c1  m notin A, m notin B, A != B
  c2  m notin A, m in B,    A != B        (possible: N=2,m=1,B={1},A={2})
  c3  A = B  (forces m notin B)           -- this is (BD) at ord_S = 0
  and the type split is orthogonal (same combinatorics, different engine).
"""
import sys, itertools
sys.path.insert(0, '.')
import q484_engine as E

def cell(A, B, m):
    if A == B:
        return 'c3-diag'
    return 'c1-mnotinB' if m not in B else 'c2-minB'

print("N  m   c1-mnotinB  c2-minB  c3-diag   TOTAL  | plant(A not>=B, m notin A)")
grand = {}
for N in range(1, 6):
    P = E.subsets(N)
    for m in range(1, N + 1):
        cnt = {'c1-mnotinB': 0, 'c2-minB': 0, 'c3-diag': 0}
        plant = 0
        for A in P:
            if m in A:
                continue
            for B in P:
                if E.geq(A, B):
                    cnt[cell(A, B, m)] += 1
                else:
                    plant += 1
        tot = sum(cnt.values())
        print(f"{N}  {m}   {cnt['c1-mnotinB']:9d}  {cnt['c2-minB']:7d}  {cnt['c3-diag']:7d}  {tot:6d}  | {plant}")
        for k, v in cnt.items():
            grand[k] = grand.get(k, 0) + v
print()
print("GRAND (N=1..5, all m, ONE type -- multiply by #eps for the full census):", grand,
      "total", sum(grand.values()))

# the N<=4 sub-domain, to reconcile with the paper's 699/699
g4 = {'c1-mnotinB': 0, 'c2-minB': 0, 'c3-diag': 0}
pl4 = 0
for N in range(1, 5):
    P = E.subsets(N)
    for m in range(1, N + 1):
        for A in P:
            if m in A:
                continue
            for B in P:
                if E.geq(A, B):
                    g4[cell(A, B, m)] += 1
                else:
                    pl4 += 1
print("N<=4 one type:", g4, "total", sum(g4.values()), " plant", pl4)
print("N<=4 x (C:2 eps) + (B:1 eps) = x3 weighted:",
      {k: v for k, v in g4.items()}, "-> C 2*%d + B 1*%d = %d"
      % (sum(g4.values()), sum(g4.values()), 3 * sum(g4.values())))
