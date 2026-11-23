"""Data corroboration of the square sandwich (Lemma sand) and its equality case, and of Lemma cls1.
For every pair of journeys in the census pickles (J_x, J_y of every diverging B/C pair and B/B pair) and every non-rhombus
M0 tile u removed by both, at steps k (by J) and j (by J'), in both orders:  e_u = Pi_T'(R'_j) - Pi_T(R_k+1).
Squares: must be <= 0, and = 0 only for the (1,0)x(1,-1) square with J' at a class-0 b-step.  Triangles: the lesion.
Also: S = H_p cap W contains the square of H_p, on every role."""
import pickle
exec(open('/home/clio/projects/proofs/code-1123-Q597/sandwich.py').read().split("bad = 0")[0])
SH2 = ('sq', ((0, 0, 0, 0), (0, 0, 1, -1), (1, 0, 0, 0), (1, 0, 1, -1)))
st = Counter()
def rh0(s): return it(next(u for u in s[1] if u[0] == 'rh'))
def rh1(s): return it(next(u for u in s[2] if u[0] == 'rh'))
for SIZE in (5, 6):
    for P in pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb')):
        M0 = set(P['M0']); Js = [r[1] for r in P['roles']]
        Ts = ['B', 'B'] if P['kind'] == 'BB' else ['B' if r[0].endswith('B') else 'C' for r in P['roles']]
        rem = [{u: k for k, s in enumerate(J) for u in s[1] if u in M0 and u[0] != 'rh'} for J in Js]
        for i in (0, 1):
            J, Jp, T, Tp = Js[i], Js[1 - i], Ts[i], Ts[1 - i]
            for u in set(rem[0]) & set(rem[1]):
                k, j = rem[i][u], rem[1 - i][u]
                e = sub(Pi(rh0(Jp[j]), Tp), Pi(rh1(J[k]), T)); uu = it(u)
                cj = (cls(Tp, rh0(Jp[j])), Jp[j][0])
                if u[0] == 'sq':
                    eq = e == (0, 0, 0, 0)
                    st[('size', SIZE, P['kind'], 'SQUARE e<=0', all(v <= 0 for v in e), 'e=0', eq,
                        'shape2 & J\' (0,b)' if eq else '', (shape(uu) == SH2 and cj == (0, 'b')) if eq else '')] += 1
                else:
                    st[('size', SIZE, P['kind'], 'TRIANGLE (lesion) e<=0', all(v <= 0 for v in e))] += 1
            # S contains the square of H_p
            W = set(rem[1 - i]); p = min(k for k, s in enumerate(J) if set(s[1]) & W)
            sq = [u for u in J[p][1] if u[0] == 'sq']
            st[('size', SIZE, P['kind'], 'square of H_p in W', len(sq) == 1 and sq[0] in W)] += 1
# Lemma cls1: every configuration removing a shape-2 square ends at a class-1 rhombus; stationary steps = class-0 b-steps
for (T, ck, typ), (Rc, rem_, Rn, cre) in CF.items():
    if any(u[0] == 'sq' and shape(u) == SH2 for u in rem_): st[('TABLE removes shape-2 square', T, cls(T, Rc), typ, 'ends class', cls(T, Rn))] += 1
    st[('TABLE stationary', T, cls(T, Rc), typ, sub(Pi(Rn, T), Pi(Rc, T)) == (0, 0, 0, 0))] += 1
for k in sorted(st, key=str): print(st[k], k)
