import sympy as sp
from sympy import symbols, expand, simplify, together
import sys
sys.path.insert(0, '.')
from extract_n5_eigenvalues import (
    syt_of_shape, build_Ti_seminormal, verify_hecke_relations
)

q = symbols('q')

lam = (2, 1)
n = 3

syts = syt_of_shape(lam)
syts_keyed = sorted(syts, key=lambda T: tuple(T[k] for k in sorted(T)))
index = {tuple((k, T[k]) for k in sorted(T)): idx
         for idx, T in enumerate(syts_keyed)}

print(f"SYTs of {lam}:")
for idx, T in enumerate(syts_keyed):
    cells = [(k, T[k]) for k in sorted(T)]
    print(f"  idx {idx}: {cells}")

mats = {}
for i in range(1, n):
    M = build_Ti_seminormal(lam, i, syts_keyed, index)
    print(f"\nT_{i} =")
    print(sp.pretty(M))
    mats[i] = M

# Verify T_i^2 = (q-1)T_i + q
for i in range(1, n):
    lhs = expand(mats[i] * mats[i])
    rhs = expand((q-1)*mats[i] + q*sp.eye(mats[i].shape[0]))
    diff = sp.Matrix([[simplify(lhs[a,b] - rhs[a,b]) for b in range(lhs.shape[1])]
                      for a in range(lhs.shape[0])])
    print(f"\nT_{i}^2 - ((q-1)T_{i} + q) =")
    print(sp.pretty(diff))

# Braid
lhs = expand(mats[1] * mats[2] * mats[1])
rhs = expand(mats[2] * mats[1] * mats[2])
diff = sp.Matrix([[simplify(lhs[a,b] - rhs[a,b]) for b in range(lhs.shape[1])]
                  for a in range(lhs.shape[0])])
print(f"\nBraid T_1 T_2 T_1 - T_2 T_1 T_2 =")
print(sp.pretty(diff))
