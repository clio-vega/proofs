"""Existence half: if J_o (toward T') stopped at RJ = R^o_{e+1}, then RJ = cell (i,j) of nest T' (Lemma grid / packed_at) and
cells (i-1,j), (i,j-1) are nest rhombi RJ-u, RJ-w when i>0, j>0.  Do Z1, Z2 (canonical) conflict with them?  Is RJ even a
T'-frame rhombus?   usage: python3 zstop.py"""
import sys
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict
FR = {'B': ((-1, 1, 0, 0), (0, 0, 1, -1)), 'C': ((1, 0, 0, 0), (0, 0, -1, 0))}
Z = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
     'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
for Tp in 'CB':
    Rs = CF[Tp][(k0, 'b')][0]; HJ, RJ, CJ = step(Tp, Rs, 'b'); o = min(Rs[1])
    zs = [('tri', frozenset(tuple(a + b for a, b in zip(P, o)) for P in z)) for z in Z[Tp]]
    u, w = FR[Tp]
    frame = any({q, tuple(a+b for a, b in zip(q, u)), tuple(a+b for a, b in zip(q, w)), tuple(a+b+c for a, b, c in zip(q, u, w))} == set(RJ[1]) for q in RJ[1])
    Nu = trs(RJ, tuple(-x for x in u)); Nw = trs(RJ, tuple(-x for x in w))
    print('J_o ->', Tp, '| RJ is T\'-frame:', frame, '| Z1 vs (i-1,j):', bool(conflict([zs[0]], [Nu])), ' Z1 vs (i,j-1):', bool(conflict([zs[0]], [Nw])),
          '| Z2 vs (i-1,j):', bool(conflict([zs[1]], [Nu])), ' Z2 vs (i,j-1):', bool(conflict([zs[1]], [Nw])))
