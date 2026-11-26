"""B/B CONTROL copy of corr.py. PROVE 11-26 c2: the corridor on M0 (Merge) from the shared rhombus R*.
Prediction (written in PREDICTION.md before running): in case b, for i >= 2 while E_i and L_i are taken at the same rhombus,
they remove the same tiles; while the common type is an a-type the steps are identical; the first b-step (index k+2) is
(0,b) for both, same removed tiles, different created.  (H2): at every shared rhombus both journeys take a step.
Also: the class-0 config table identity (a-types identical for B and C, b-type same removed set).
usage: python3 corr.py SIZE"""
import sys, glob, pickle
from collections import Counter
SIZE = int(sys.argv[1])
def rh(S): return next(t for t in S if t[0] == 'rh')
st = Counter()
for fn in sorted(glob.glob('/tmp/step0BB_%d_*.pkl' % SIZE)):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        if case != 'b':
            # case a: E_2 and L_2 at different rhombi?
            same = len(AE) > pE+2 and len(AL) > pL+2 and rh(AE[pE+2][1]) == rh(AL[pL+2][1])
            st[('a', 'E_2, L_2 at the same rhombus', same)] += 1
            continue
        i = 2; kk = None
        while True:
            hasE, hasL = len(AE) > pE+i, len(AL) > pL+i
            if not (hasE and hasL):
                st[('b', 'H2 violated: a journey stops at a shared rhombus', 'i', i, hasE, hasL)] += 1; break
            e, l = AE[pE+i], AL[pL+i]
            if rh(e[1]) != rh(l[1]): st[('b', 'rhombi differ at i', i)] += 1; break
            same_rem = set(e[1]) == set(l[1]); ident = same_rem and set(e[2]) == set(l[2])
            if e[0] == 'b' or l[0] == 'b':
                st[('b', 'first b at shared rhombus: types', e[0], l[0], 'same removed', same_rem, 'identical', ident)] += 1
                kk = i - 2; break
            st[('b', 'shared a-step: same removed', same_rem, 'identical', ident)] += 1
            i += 1
        if kk is not None:
            st[('b', 'k (corridor length) == m_L - 3', kk == laL - 3)] += 1
            # after the exit: does any later step of J_E equal a later step of J_L?
            later = {(s[0], s[1], s[2]) for s in AL[pL+kk+3:]}
            st[('b', 'a later E step equals a later L step', any((s[0], s[1], s[2]) in later for s in AE[pE+kk+3:]))] += 1
print('SIZE', SIZE)
for k in sorted(st, key=str): print(st[k], k)
