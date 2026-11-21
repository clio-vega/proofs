"""PROVE 11-21 c2: controls for shpast.py.
(C1) same target (T = T'): S_J for T-words ending in the class-0 b at R*, S_A for the T-continuation. Since the same-target
     journeys shadow on through the b-step, this is the B/B situation; does the detector fire?
(C2) S_A translated by each of the 12 unit vectors of Lambda_0 + Lambda_1: how many translates meet S_J? (adjacency, liveness)
(F)  the FUTURE case (J removes t after e): feasible tiles with psi_T <= , lam >= (A) and psi_T' <= , lam >= (J) in boxes
     of growing size: if the count grows with the box and hits the boundary, no window bound exists.
usage: python3 shpast_ctl.py"""
import sys
sys.argv = ['shpast.py', '6']
exec(open('/home/clio/projects/proofs/code-1121c2-Q597/shpast.py').read().split("if __name__ == '__main__':")[0])
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W12 = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
for T in 'BC':
    Sj, nj = SJ(T, 8); Sa, na, split = SA(T, 5)
    print('C1 same target', T, '|S_J|', len(Sj), '|S_A|', len(Sa), 'S_J & S_A', len(Sj & Sa))
for T, Tp in (('B', 'C'), ('C', 'B')):
    Sj, nj = SJ(Tp, 8); Sa, na, split = SA(T, 5)
    hits = sum(1 for v in W12 if {trs(t, v) for t in Sa} & Sj)
    print('C2', T, 'unit translates of S_A meeting S_J:', hits, 'of 12')
def future_feasible(T, Tp, box):
    n = ob = 0
    for s in shapes:
        for v in itertools.product(range(-box, box + 1), repeat=4):
            t = trs(s, v); oA, oJ = offs(t, Rstar, T), offs(t, Rstar, Tp)
            okA = any(oA[0] <= q[0] and oA[1] <= q[1] and oA[2] >= q[2] and oA[3] >= q[3] for q in QS[T])
            okJ = any(oJ[0] <= q[0] and oJ[1] <= q[1] and oJ[2] >= q[2] and oJ[3] >= q[3] for q in QS[Tp])
            if okA and okJ:
                n += 1; ob += max(abs(x) for x in v) == box
    return n, ob
for box in (3, 5, 7):
    print('F future feasible (B,C) box', box, future_feasible('B', 'C', box))
