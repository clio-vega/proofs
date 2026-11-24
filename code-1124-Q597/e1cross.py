"""Cross test at the split (abstract configurations): does J_o's step-e+1 configuration of each type at RJ conflict with
A's created tiles CA (+ R#)?  And with J_o's own CJ (sanity: must be 0 for every placeable type)?
usage: python3 e1cross.py"""
import sys
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict
nr = lambda S: frozenset(u for u in S if u[0] != 'rh')
for T in 'BC':
    Tp = OTHER[T]; Rs = CF[T][(k0, 'b')][0]
    HA, Rsh, CA = step(T, Rs, 'b'); HJ, RJ, CJ = step(Tp, Rs, 'b')
    for tJ in ('a0', 'a1', 'b'):
        H1, _, _ = step(Tp, RJ, tJ)
        print('T', T, 'J_o(->%s) e+1 type (1,%s)' % (Tp, tJ), '| conflicts with CJ+RJ:', len(conflict(H1 - {RJ}, CJ)),
              '| with A side CA+R#:', [(a[0], b[0]) for a, b in conflict(H1, CA | {Rsh})], '| inside old hexagon:', len(conflict(H1 - {RJ}, HJ)))
