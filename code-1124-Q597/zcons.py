"""Lemma Zblock(1): Z1 (and Z2) meet no tile of C_J + R^o_{e+1}, of C_A + R#, or of the split hexagon H* in interior
(canonical configuration).  usage: python3 zcons.py"""
import sys
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict
Z = {'C': [((0, 0, 0, -1), (0, 0, 0, 0), (0, 0, 1, -1)), ((0, 0, 0, 0), (0, 0, 1, -1), (0, 0, 1, 0))],
     'B': [((0, 1, 0, 1), (0, 2, 0, 1), (1, 1, 0, 1)), ((0, 1, 0, 1), (1, 0, 0, 1), (1, 1, 0, 1))]}
for Tp in 'CB':
    T = OTHER[Tp]; Rs = CF[Tp][(k0, 'b')][0]; HJ, RJ, CJ = step(Tp, Rs, 'b'); HA, Rsh, CA = step(T, Rs, 'b'); o = min(Rs[1])
    for n, z in enumerate(Z[Tp]):
        zt = ('tri', frozenset(tuple(a + b for a, b in zip(P, o)) for P in z))
        print('J_o->', Tp, 'Z%d' % (n + 1), '| vs CJ+RJ:', len(conflict([zt], CJ | {RJ})), '| vs CA+R#:', len(conflict([zt], CA | {Rsh})),
              '| vs H*:', len(conflict([zt], HJ)), '| in H*:', zt in HJ)
