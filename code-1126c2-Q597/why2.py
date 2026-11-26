"""Name the killing constraint of local2.py, one role at a time; match the survivors to the census classes (joint5.pkl,
joint6.pkl); B/B control (T_E = T_L).
For each contact pair and each candidate type t of E_1 alone: is G_t(R_E^1) compatible with the M0 tiles of the two
contact hexagons?  If exactly one t survives for each role, (Share_2) is a statement about two forced tile sets."""
import sys, itertools, pickle, os
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
aT = {'B': 'a1', 'C': 'a0'}
def sqof(S): return next(t for t in S if t[0] == 'sq')
def shp(t): b = base(t); return (t[0], tuple(sorted(tuple(a - c for a, c in zip(P, b)) for P in t[1])))
def configs_at(T, R): return {typ: step(T, R, typ) for typ in TYPES if (T, classkey(R), typ) in CF}
def killers(S, U):
    """pairs (s in S, u in U) that overlap without being equal"""
    return [(s[0], u[0]) for s in S for u in U if s != u and overlap(s, u)]
def contacts(TL, TE, Lcls=1, Ltyp=None):
    Ltyp = Ltyp or aT[TL]
    RcL, remL, RnL, creL = CF[(TL, K1[TL] if Lcls == 1 else K0, Ltyp)]
    for (T_, ckE, tE0), (RE, remE, RnE, creE) in CF.items():
        if T_ != TE: continue
        if cls(TE, RE) == 1 and tE0 == aT[TE]: continue
        if TE == TL and (ckE, tE0) == (classkey(RcL), Ltyp): continue
        u = sqof(remE); uL = sqof(remL)
        if shp(u) != shp(uL): continue
        D = tuple(a - b for a, b in zip(base(u), base(uL)))
        yield (cls(TE, RE), tE0), RE, set(remE) | {RE}, set(creE) | {RnE}, RnE, \
              tr(RcL, D), {tr(t, D) for t in remL} | {tr(RcL, D)}, {tr(t, D) for t in creL} | {tr(RnL, D)}, tr(RnL, D), u
MODE = sys.argv[1] if len(sys.argv) > 1 else 'BC'
REPS = []
for TL in 'BC':
    TE = OTHER[TL] if MODE == 'BC' else TL
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        nHE, nHL = set(nonrh(HE)), set(nonrh(HL))
        print('CONTACT T_L', TL, 'T_E', TE, 'E0', cE0, 'hexagons compatible', not killers(nHE, nHL))
        cfE, cfL = configs_at(TE, RE1), configs_at(TL, RL1)
        okE, okL = [], []
        for t, (rem, _, cre) in cfE.items():
            G = set(nonrh(rem)); k = killers(G, nHE | nHL)
            print('   E1 type', t, 'G vs M0 hexagon tiles: overlapping pairs', k or 'none', '| G meets nH_E as equal tile', bool(G & nHE))
            if not k and not (G & nHE): okE.append(t)
        for t, (rem, _, cre) in cfL.items():
            G = set(nonrh(rem)); k = killers(G, nHE | nHL)
            print('   L1 type', t, 'G vs M0 hexagon tiles: overlapping pairs', k or 'none', '| G meets nH_L as equal tile', bool(G & nHL))
            if not k and not (G & nHL): okL.append(t)
        print('   forced: E1', okE, 'L1', okL)
        for tE, tL in itertools.product(okE, okL):
            GE, GL = set(nonrh(cfE[tE][0])), set(nonrh(cfL[tL][0]))
            OE, OL = GE - nHE - nHL, GL - nHE - nHL
            print('   PAIR', tE, tL, 'GE vs GL overlapping pairs', killers(GE, GL) or 'none', 'Oth(E1)==Oth(L1)', OE == OL,
                  'Oth kinds', sorted(t[0] for t in OE), sorted(t[0] for t in OL))
            if killers(GE, GL): continue
            b = min(u[1]); f = lambda S: frozenset((t[0], frozenset(tuple(a - c for a, c in zip(P, b)) for P in t[1])) for t in S)
            REPS.append((TL, (f(HE), f(CE0), f(set(cfE[tE][0])), f(set(cfE[tE][2])), f(HL), f(CL0), f(set(cfL[tL][0])), f(set(cfL[tL][2])))))
if MODE == 'BC':
    for SIZE in (5, 6):
        fn = '/home/clio/projects/proofs/code-1126c2-Q597/joint%d.pkl' % SIZE
        if not os.path.exists(fn): print('no', fn); continue
        J = pickle.load(open(fn, 'rb'))
        local = {r for r in REPS}
        hit = sum(n for (case, TL, tup), n in J.items() if (TL, tup) in local)
        print('CENSUS size', SIZE, 'records', sum(J.values()), 'classes', len(J), 'records in a local survivor class', hit,
              'census classes not local', sum(1 for (c, TL, tup) in J if (TL, tup) not in local),
              'local classes not seen', len(local - {(TL, tup) for (c, TL, tup) in J}))
