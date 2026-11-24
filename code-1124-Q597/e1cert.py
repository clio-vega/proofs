"""PROVE 11-24 main finite check for (E1).  Abstract configurations (words.configs via the hsderive preamble, as compi.py).
R* = canonical class-0 rhombus; A (toward T) and J_o (toward T') both take (0,b) at R* with the SAME hexagon (asserted).
R# = A's image, RJ = J_o's image R^o_{e+1}; CA, CJ = created tiles.
(N) non-last: A takes a class-1 step t_A at R#.  Its non-rh tiles lie in M0 & Mx, so they are on J_o's board at e+1.
    For each t_A: is nr(H_A(t_A)) compatible with J_o's tiles CJ (else t_A is impossible)?  For each t_J at RJ:
    does the t_J hexagon conflict (interior overlap, distinct) with nr(H_A(t_A))?   Claim: conflict <=> t_J = straight.
(L) last: A stops at R#.  Then R# = cell (i,j) of nest T and cells (i-1,j), (i,j-1) are nest rhombi R#-u_T, R#-w_T
    when i>0, j>0.  Does the t_J hexagon conflict with them?
usage: python3 e1cert.py"""
import sys
src = open('/home/clio/projects/proofs/code-1121c2-Q597/hsderive.py').read().split('\nCCS = ')[0]
exec(src)
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import conflict, overlap
STR = {'B': 'a1', 'C': 'a0'}
FR = {'B': ((-1, 1, 0, 0), (0, 0, 1, -1)), 'C': ((1, 0, 0, 0), (0, 0, -1, 0))}
nr = lambda S: frozenset(u for u in S if u[0] != 'rh')
def neg(v): return tuple(-x for x in v)
def trystep(T, R, t):
    try: return step(T, R, t)
    except (KeyError, AssertionError): return None
for T in 'BC':
    Tp = OTHER[T]
    Rs = CF[T][(k0, 'b')][0]; assert Rs == CF[Tp][(k0, 'b')][0]
    HA, Rsh, CA = step(T, Rs, 'b'); HJ, RJ, CJ = step(Tp, Rs, 'b'); assert HA == HJ
    print('==== T =', T, ' J_o toward', Tp, ' straight of J_o = (1,%s)' % STR[Tp])
    print('  R# == RJ ?', Rsh == RJ, '| CA == CJ ?', CA == CJ, '| conflict(CA,CJ) pairs:', len(conflict(CA, CJ)))
    Jbd = CJ | {RJ}
    for tA in ('a0', 'a1', 'b'):
        s = trystep(T, Rsh, tA)
        if s is None: print('  (N) t_A', tA, 'not placeable'); continue
        TA = nr(s[0])
        compatJ = not conflict(TA, Jbd); selfA = not conflict(TA, CA - TA)
        line = []
        for tJ in ('a0', 'a1', 'b'):
            sj = trystep(Tp, RJ, tJ)
            if sj is None: line.append((tJ, 'np')); continue
            line.append((tJ, len(conflict(sj[0], TA))))
        print('  (N) t_A=(1,%s): A-tiles in CA: %d | compatible with J_o tiles at e+1: %s | J_o e+1 conflicts:' % (tA, len(TA & CA), compatJ), line)
    u, w = FR[T]
    # check R# is a T-frame rhombus
    P0 = min(Rsh[1]); frame = any({q, tuple(a+b for a, b in zip(q, u)), tuple(a+b for a, b in zip(q, w)), tuple(a+b+c for a, b, c in zip(q, u, w))} == set(Rsh[1]) for q in Rsh[1])
    Nu = trs(Rsh, neg(u)); Nw = trs(Rsh, neg(w))
    print('  (L) R# is a T-frame rhombus:', frame, '| nest neighbours compatible with CA:', not conflict({Nu, Nw}, CA | {Rsh}), 'with J_o tiles CJ+RJ:', not conflict({Nu, Nw}, Jbd))
    for tJ in ('a0', 'a1', 'b'):
        sj = trystep(Tp, RJ, tJ)
        if sj is None: print('    t_J', tJ, 'np'); continue
        print('    (L) t_J=(1,%s): conflicts with cell (i-1,j) rhombus R#-u: %d, with (i,j-1) rhombus R#-w: %d' % (tJ, len(conflict(sj[0], {Nu})), len(conflict(sj[0], {Nw}))))
