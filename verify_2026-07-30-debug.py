"""Debug: what is M_{(3,0,0)} really, and does it match v * P?"""
import sys
sys.path.insert(0, '/home/clio/projects/scratch/2026-07-24-pieri-affine-CR')
sys.path.insert(0, '/home/clio/projects/proofs')

import sympy as sp
from sympy import expand, simplify
from verify_2026_07_30_interior_identity import (
    build_M_iterated, is_interior, is_dominant, in_Pc, V_coeff
)
from hall_littlewood_engine import hall_littlewood_P, t_sym, v_lambda
from demazure_engine import make_vars

n, c, max_size = 3, 3, 4
X = make_vars(n)
M_dict, status_dict, P_dict, v_dict = build_M_iterated(n, c, max_size)

print("M_{(3,0,0)}:")
print("  computed by script:", sp.expand(M_dict[(3,0,0)]))
print("  v * P at (3,0,0):  ", sp.expand(v_dict[(3,0,0)] * P_dict[(3,0,0)]))
print("  difference:         ", sp.simplify(M_dict[(3,0,0)] - v_dict[(3,0,0)] * P_dict[(3,0,0)]))

print("\nM_{(3,1,0)}:")
print("  computed:", sp.expand(M_dict[(3,1,0)]))
print("  v * P:   ", sp.expand(v_dict[(3,1,0)] * P_dict[(3,1,0)]))
print("  diff:    ", sp.simplify(M_dict[(3,1,0)] - v_dict[(3,1,0)] * P_dict[(3,1,0)]))

# Compare to what I derived from Pieri (2,0,0) manually
print("\nHand derivation of M_{(3,0,0)} from Pieri (2,0,0):")
print("  P_{(3,0,0)} - t * P_{(2,1,0)} =", sp.expand(P_dict[(3,0,0)] - t_sym * P_dict[(2,1,0)]))
