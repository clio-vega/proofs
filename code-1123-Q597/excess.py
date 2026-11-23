"""Per-role excess e(s) = Pi_T'(R^o_c0) - Pi_T(R_p+1) on (i) the Theorem-local survivors (B/C: ccsurv.pkl; B/B control:
ccsurv-same.pkl, T' = T) and (ii) every diverging role of the census pickles (size 5 full, size 6 KEEP=3), computed from
the real journeys.  Also the Share consequence e'(x)+e'(y) = 0 on real pairs, and (Share) itself.  usage: python3 excess.py"""
import os, pickle
exec(open('/home/clio/projects/proofs/code-1123-Q597/commonpot.py').read().split("st = Counter()")[0])
st = Counter()
for tag in ('', '-same'):
    FF = pickle.load(open('/home/clio/projects/proofs/code-1123-Q597/ccsurv%s.pkl' % tag, 'rb'))
    for s in FF:
        T = s['T']; To = T if tag else OTHER[T]
        RA = s['RA']; R1 = step(T, RA, s['t'])[1]; Rc = s['seg'][s['c0']][0]
        e = sub(Pi(Rc, To), Pi(R1, T))
        st[('LOCAL', 'B/B' if tag else 'B/C', T, 'e<=0', all(v <= 0 for v in e), 'e!=0', any(e), 'e', e)] += 1
def rh0(s): return it(next(u for u in s[1] if u[0] == 'rh'))
def rh1(s): return it(next(u for u in s[2] if u[0] == 'rh'))
for SIZE in (5, 6):
    PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
    for P in PAIRS:
        M0 = set(P['M0']); E = []; EP = []
        Js = [r[1] for r in P['roles']]
        Ts = ['B', 'B'] if P['kind'] == 'BB' else ['B' if r[0].endswith('B') else 'C' for r in P['roles']]
        rem = [{u: k for k, s in enumerate(J) for u in s[1] if u in M0} for J in Js]
        for i in (0, 1):
            J, Jo, T, To = Js[i], Js[1 - i], Ts[i], Ts[1 - i]
            W = set(rem[1 - i])
            p = min(k for k, s in enumerate(J) if set(s[1]) & W)            # first contact (= divergence, Thm fc)
            Hp = set(J[p][1]); c0 = min(rem[1 - i][u] for u in Hp & W)
            e = sub(Pi(rh0(Jo[c0]), To), Pi(rh1(J[p]), T)); ep = sub(Pi(rh0(Jo[c0]), To), Pi(rh0(J[p]), T))
            E.append(e); EP.append(ep)
            st[('DATA size', SIZE, P['kind'], 'e<=0', all(v <= 0 for v in e), 'e!=0', any(e))] += 1
            pj = min(k for k, s in enumerate(Jo) if set(s[1]) & set(rem[i]))
            st[('DATA size', SIZE, P['kind'], 'c0 - p_other', c0 - pj)] += 1
        st[('DATA size', SIZE, P['kind'], "e'(x)+e'(y)", add(EP[0], EP[1]), 'e(x)+e(y)', add(E[0], E[1]))] += 1
for k in sorted(st, key=str): print(st[k], k)
