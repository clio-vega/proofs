"""Case a, the same tools as case b's back window.  Proved inputs: E_0, E_1=(0,b), L_0=(1,a_TL), L_1=(1,a_TE) (Thm share2 + h1),
E^L_0, E^L_1=(1,b) (Prop rerun), L^E_0 = L_1-alpha, L^E_1 = L_0+gamma (Thm bsall, C=(L_1)).
Checks: (N) no stop of J_E at R^E_2: an L^E step removes a tile overlapping R^E_2.  (P) pair lemma on Y before E^L_2:
present = P1 (G(E_2) minus tiles removed by L_0,L_1,E^L_0,E^L_1) u Kcur (K updated by E^L_0,E^L_1) u G(E^L_2).
(R) resync after step 2.  (S) stale M0 tiles: G(L_0)uG(L_1) and G(E^L_2), not removed by E_0..E_2: offset past Q at R^E_2/R*?
(W) the middles (E_0..2, E^L_0..2, L_0..1, L^E_0..1) against windows.pkl W['a']."""
import itertools, pickle
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
exec(open('/home/clio/projects/proofs/code-1127-Q597/r1.py').read().split("H = Counter()")[0].split('exec(')[0])
ETA = {'B': lambda a, b: b, 'C': lambda a, b: a + b}
def psi(P, T): return (min(-ETA[T](p[2], p[3]) for p in P), min(ETA[T](p[0], p[1]) for p in P))
def Af(P): return (min(p[0] for p in P), min(p[2] for p in P))
def om(t, R, T):
    a, b = psi(t[1], T), psi(R[1], T); c, d = Af(t[1]), Af(R[1]); return (a[0]-b[0], a[1]-b[1], c[0]-d[0], c[1]-d[1])
Q = {}
for T in 'BC':
    for c in (0, 1):
        R = canon(T, c)
        for typ, (rem, Rn, cre) in configs_at(T, R).items():
            for g in nonrh(rem): Q.setdefault(T, set()).add(om(g, R, T))
def past(f, T):
    q = Q[T]
    return f[0] > max(x[0] for x in q) or f[1] > max(x[1] for x in q) or f[2] < min(x[2] for x in q) or f[3] < min(x[3] for x in q)
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
def pairwise_ok(S): return all(a == b or not overlap(a, b) for a, b in itertools.combinations(list(S), 2))
WIN = pickle.load(open('/home/clio/projects/proofs/code-1126-Q597/windows.pkl', 'rb'))
REPS = {key: v for key, v in WIN["a"].items()}
def conv(typ, R, res):
    rem, Rn, cre = res
    F = lambda S: frozenset((t[0], frozenset(t[1])) for t in S)
    return (typ, F(set(rem) | {R}), F(set(cre) | {Rn}))
def shift(s, o):
    f = lambda S: frozenset((t[0], frozenset(tuple(a - b for a, b in zip(p, o)) for p in t[1])) for t in S)
    return (s[0], f(s[1]), f(s[2]))
matched = set(); ncls = 0
for TL in 'BC':
    TE = OTHER[TL]; al = ALPHA[TL]; mal = tuple(-a for a in al)
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        cfE, cfL = configs_at(TE, RE1), configs_at(TL, RL1)
        if 'b' not in cfE or aT[TE] not in cfL: continue
        remE1, RE2, creE1 = cfE['b']; remL1, RL2, creL1 = cfL[aT[TE]]
        nHE, nHL = set(nonrh(HE)), set(nonrh(HL)); GE1, GL1 = set(nonrh(remE1)), set(nonrh(remL1))
        if not pairwise_ok(nHE | nHL | GE1 | GL1): continue      # the (Share2) classification, case a
        ncls += 1
        CL1 = set(creL1)
        KK = set(nonrh(CL0 | CL1)) | (nHE - nHL - GL1)
        c0 = configs_at(TE, RE); pres = [t for t in c0 if set(nonrh(c0[t][0])) <= KK]
        assert len(pres) == 1
        EL0 = c0[pres[0]]; rem0, R1, cre0 = EL0
        KA = (KK - set(nonrh(rem0))) | set(nonrh(cre0))
        c1 = configs_at(TE, R1); pres1 = [t for t in c1 if set(nonrh(c1[t][0])) <= KA]; assert pres1 == ['b']
        EL1 = c1['b']; rem1, REL2, cre1 = EL1
        KB = (KA - set(nonrh(rem1))) | set(nonrh(cre1))
        # L^E: L^E_0 = L_1 - alpha at RL (common prefix); L^E_1 = L_0's type at the rhombus L^E_0 creates
        assert tr(RL1, mal) == RL, 'L^E_0 at L contact rhombus'
        LE0 = step(TL, RL, aT[TE]); assert set(LE0[0]) == {tr(g, mal) for g in remL1}
        LE1 = step(TL, LE0[1], aT[TL])
        line = dict(TL=TL, E0=cE0, EL0=pres[0], RE2cls=cls(TE, RE2), REL2cls=cls(TE, REL2))
        # (N)
        line['N: L^E_0 tile over R^E_2'] = [g[0] for g in nonrh(LE0[0]) if overlap(g, RE2)]
        line['N: L^E_1 tile over R^E_2'] = [g[0] for g in nonrh(LE1[0]) if overlap(g, RE2)]
        line['N: R^E_2 != L^E rhombi'] = RE2 not in (RL, LE0[1])
        # (P)
        GL01 = nHL | GL1
        res = {}
        for x, rx in configs_at(TE, RE2).items():
            P1 = set(nonrh(rx[0])) - GL01
            for y, ry in configs_at(TE, REL2).items():
                ok = pairwise_ok(P1 | KB | set(nonrh(ry[0])))
                if ok: res[(x, y)] = (ry[1] == rx[1], ry, rx)
            line.setdefault('G(E_2) cap G(L0 u L1) by x', {})[x] = len(set(nonrh(rx[0])) & GL01)
        line['P: compatible pairs'] = sorted(res)
        line['R: same rhombus after'] = [res[p][0] for p in sorted(res)]
        # (S) stale M0: G(L0)uG(L1) not removed by E_0,E_1,E_2 ; G(E^L_2) M0 tiles (not in KB's created part) not removed by E_0..E_2
        st = []
        for (x, y), (_, ry, rx) in sorted(res.items()):
            GE012 = set(nonrh(HE)) | GE1 | set(nonrh(rx[0]))
            Rs = rx[1]
            cand = [('L', g) for g in GL01 - GE012] + [('EL2', g) for g in set(nonrh(ry[0])) - GE012 - set(nonrh(CL0 | CL1)) - set(nonrh(cre0)) - set(nonrh(cre1))]
            st.append((x, y, [(w, g[0], past(om(g, Rs, TE), TE)) for w, g in cand]))
            # (W)
            mids = (tuple(conv(t[0], R_, r_) for t, R_, r_ in ((cE0[1], RE, c0[cE0[1]]), ('b', RE1, cfE['b']), (x, RE2, rx))),
                    tuple(conv(t[0], R_, r_) for t, R_, r_ in ((pres[0], RE, EL0), ('b', R1, EL1), (y, REL2, ry))),
                    tuple(conv(t[0], R_, r_) for t, R_, r_ in ((aT[TL], RL, step(TL, RL, aT[TL])), (aT[TE], RL1, cfL[aT[TE]]))),
                    tuple(conv(t[0], R_, r_) for t, R_, r_ in ((aT[TE], RL, LE0), (aT[TL], LE0[1], LE1))))
            sq = [t for t in HE if t[0] == 'sq' and t in HL]
            o = min(sq[0][1]) if sq else None
            if o is not None:
                key = (TL, tuple(tuple(shift(s, o) for s in m) for m in mids))
                hit = [v[1] for k_, v in REPS.items() if k_ == key]
                if hit: matched.add(key)
                st[-1] = st[-1] + (('W match', hit),)
        line['S: stale (who, kind, past at R*)'] = st
        print(line)
print('case-a classes', ncls, '| census reps', len(REPS), '| matched', len(matched))
