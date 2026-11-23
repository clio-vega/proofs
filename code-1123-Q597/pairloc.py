"""PROVE 11-23, Q597 (Share) as a finite check at the PAIR level.
Failure of (Share) (Lemma hs): c0(x) > p_y and c0(y) > p_x, H^x cap H^y = {} (a shared tile would be an M0 tile removed at
p_x and p_y).  Role x is a survivor sx of Theorem local for T=B (A=J_x, J_o=J_y; its segment = J_y's steps from c0(x));
role y a survivor sy for T=C translated by an unknown D in Z^4 (its segment = J_x's steps from c0(y)).
Necessary conditions (tags; proofs in the tex):
  BOX  (MU and lambda-monotone, strict because c0(x) >= p_y+1): mu_C(R^y_{c0x}) <= mu_C(R^y_{py+1}),
       lam(R^y_{c0x}) >= lam(R^y_{py+1}); same for J_x with B.  Linear in D, so D ranges over an explicit finite box.
  DIS  H^x, H^y share no tile and are interior-disjoint (both in the tiling Z = M0.J_x[:p_x].J_y[:p_y]).
  SIG  (swept regions) every tile known on X^x_{p_x} and not in W_y -- (H^x minus S_x) and comG_x -- is interior-disjoint
       from, and unequal to, every tile of every J_y hexagon we know (H^y step and seg_x, removed and created). Sym.
  BRD  J_y's board at p_y holds H^y and every M0 tile J_y removes later (Wseg_x); J_y's created non-rh tiles are final
       (no revisit) and are compatible with everything J_y holds later; no re-creation. Sym for J_x.
  UQX  Lemma U on J_y's board at p_y (known: H^y, Wseg_x): no other good C-config at R^y fully known. Sym.
  UQY  Lemma U on the re-run board Y^x_{p_x} = M0.J_y.J_x[:p_x], which holds R^x, G_{t2x}, H^x minus S_x and every
       non-rh tile J_y creates (final, outside J_x's prefix region): compatible, and no config other than t2x. Sym.
MODE=fail (default) or MODE=share (control: D is forced by seg_x[0] == y's contact step and seg_y[0] == x's).
LESION=<tag> drops one tag.  usage: python3 pairloc.py"""
import sys, os, itertools, pickle
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
MODE = os.environ.get('MODE', 'fail'); LES = os.environ.get('LESION', '')
F = pickle.load(open('/home/clio/projects/proofs/code-1123-Q597/ccsurv.pkl', 'rb'))
SB = [s for s in F if s['T'] == 'B']; SC = [s for s in F if s['T'] == 'C']
ETA = {'B': lambda a, b: b, 'C': lambda a, b: a + b}
def mu(t, T):
    e = ETA[T]; return (min(-e(P[2], P[3]) for P in t[1]), min(e(P[0], P[1]) for P in t[1]))
def lam(t): return (min(P[0] for P in t[1]), min(P[2] for P in t[1]))
def shmu(D, T): e = ETA[T]; return (-e(D[2], D[3]), e(D[0], D[1]))
# check the shift formulas against direct evaluation
import random; random.seed(1)
for _ in range(200):
    s = random.choice(F); t = random.choice(s['seg'][0][2]); D = tuple(random.randint(-5, 5) for _ in range(4))
    for T in 'BC':
        assert mu(tr(t, D), T) == tuple(a + b for a, b in zip(mu(t, T), shmu(D, T)))
        assert offs(tr(t, D), s['RA'], T)[:2] == tuple(a - b for a, b in zip(mu(tr(t, D), T), mu(s['RA'], T)))
    assert lam(tr(t, D)) == (lam(t)[0] + D[0], lam(t)[1] + D[2])
st = Counter(); SURV = []
def ok_pairs(S):
    return all(a == b or not overlap(a, b) for a, b in itertools.combinations(list(S), 2))
def disjoint_from(K, J):
    return all(k != j and not overlap(k, j) for k in K for j in J)
def role(s, T, D):
    RA = tr(s['RA'], D); rem, R1, cre = step(T, RA, s['t'])
    seg = [(tr(R, D), typ, [tr(x, D) for x in r], [tr(x, D) for x in c]) for (R, typ, r, c) in s['seg']]
    S = {tr(x, D) for x in s['S']}; V = {tr(x, D) for x in s['V']}
    G2 = nonrh(step(T, RA, s['t2'])[0]); comG = [g for g in G2 if g not in V]
    return dict(RA=RA, rem=rem, R1=R1, cre=cre, seg=seg, S=S, V=V, G2=G2, comG=comG, t=s['t'], t2=s['t2'],
                Wseg={x for st_ in seg for x in nonrh(st_[2])}, Cseg={x for st_ in seg for x in nonrh(st_[3])},
                segall=[x for st_ in seg for x in st_[2] + st_[3]])
def checks(x, y, Tx, Ty):
    """x: a role whose J_o is y's journey; returns the failing tag or None. One direction; called both ways."""
    Kx = [t for t in x['rem'] if t not in x['S']] + x['comG']
    Jy = y['rem'] + y['cre'] + x['segall']          # tiles of J_y hexagons known: y's contact step and x's segment
    if 'SIG' != LES and not disjoint_from(Kx, Jy): return 'SIG'
    if MODE == 'fail' and 'BRD' != LES:
        # J_y's board at p_y: H^y and the later-removed M0 tiles Wseg_x
        if x['Wseg'] & set(y['rem']): return 'BRD'
        if not ok_pairs(set(y['rem']) | x['Wseg']): return 'BRD'
        ycre = set(nonrh(y['cre']))
        if ycre & x['Wseg']: return 'BRD'                       # no revisit
        if x['Cseg'] & set(nonrh(y['rem'])): return 'BRD'      # no re-creation
        for st_ in x['seg']:                                    # final created tiles of J_y vs everything J_y holds later
            if not all(a == b or not overlap(a, b) for a in ycre for b in st_[2] + st_[3]): return 'BRD'
    if MODE == 'fail' and 'UQX' != LES:
        known = set(y['rem']) | x['Wseg']
        for o in TYPES:
            if o != y['t'] and all(g in known for g in nonrh(step(Ty, y['RA'], o)[0])): return 'UQX'
    if 'UQY' != LES:
        Yk = {x['RA']} | set(x['G2']) | {t for t in x['rem'] if t not in x['S']} | set(nonrh(y['cre'])) | x['Cseg']
        if not ok_pairs(Yk): return 'UQY'
        for o in TYPES:
            if o != x['t2'] and all(g in Yk for g in nonrh(step(Tx, x['RA'], o)[0])): return 'UQY'
    return None
Z = (0, 0, 0, 0)
for sx in SB:
    for sy in SC:
        x = role(sx, 'B', Z); y0 = role(sy, 'C', Z)
        if MODE == 'share':
            # control: c0(x) = p_y and c0(y) = p_x, i.e. seg_x[0] is y's contact step and seg_y[0] is x's
            D = tuple(a - b for a, b in zip(base(x['seg'][0][0]), base(y0['RA'])))
            Ds = [D] if tr(y0['RA'], D) == x['seg'][0][0] and x['seg'][0][1] == sy['t'] else []
        else:
            # BOX: rhombus monotonicity, strict.  J_y: R^y_{c0x} = x.seg[0].R (fixed), R^y_{py+1} = y0.R1 + D.
            #      J_x: R^x_{c0y} = y0.seg[0].R + D, R^x_{px+1} = x.R1 (fixed).
            Ryc, Ry1, Rxc, Rx1 = x['seg'][0][0], y0['R1'], y0['seg'][0][0], x['R1']
            # lam: lam(Ryc) >= lam(Ry1)+(D0,D2) ; lam(Rxc)+(D0,D2) >= lam(Rx1)
            D0hi, D2hi = lam(Ryc)[0] - lam(Ry1)[0], lam(Ryc)[1] - lam(Ry1)[1]
            D0lo, D2lo = lam(Rx1)[0] - lam(Rxc)[0], lam(Rx1)[1] - lam(Rxc)[1]
            Ds = []
            for D0 in range(D0lo, D0hi + 1):
                for D2 in range(D2lo, D2hi + 1):
                    for D1 in range(-30, 31):
                        for D3 in range(-30, 31):
                            D = (D0, D1, D2, D3)
                            mC = tuple(a + b for a, b in zip(mu(Ry1, 'C'), shmu(D, 'C')))
                            mB = tuple(a + b for a, b in zip(mu(Rxc, 'B'), shmu(D, 'B')))
                            if all(a <= b for a, b in zip(mu(Ryc, 'C'), mC)) and all(a <= b for a, b in zip(mB, mu(Rx1, 'B'))):
                                assert abs(D1) < 30 and abs(D3) < 30, 'box not bounded'
                                Ds.append(D)
            st[('box size', len(Ds))] += 1
        for D in Ds:
            y = role(sy, 'C', D)
            st['candidates'] += 1
            if MODE == 'fail' and 'DIS' != LES:
                if set(x['rem']) & set(y['rem']) or not ok_pairs(set(x['rem']) | set(y['rem'])): st['prune DIS'] += 1; continue
            f = checks(x, y, 'B', 'C') or checks(y, x, 'C', 'B')
            if f: st['prune ' + f] += 1; continue
            key = ('x', sx['cA'], sx['t'], sx['t2'], tuple((cls('C', s[0]), s[1]) for s in sx['seg']),
                   'y', sy['cA'], sy['t'], sy['t2'], tuple((cls('B', s[0]), s[1]) for s in sy['seg']), 'D', D,
                   'shared', tuple(sorted(t[0] for t in set(x['rem']) & set(y['rem']))))
            st[('SURV',) + key] += 1
print('MODE', MODE, 'LESION', LES or 'none')
for k in sorted(st, key=str): print(st[k], k)
