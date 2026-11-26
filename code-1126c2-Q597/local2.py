"""PROVE 2026-11-26 c2: the local classification behind (Share_2).
Inputs, all proved upstream:  Thm share (H_E, H_L have the same square u, an M0 tile);  Prop orient (L's contact step is
(1, a_{T_L}) of T_L, E's is not (1, a_{T_E}));  Lemma word (every step is a translate of a good configuration);  the
boards are tilings.  Hypothesis (H1): J_E and J_L both take a step at p+1 (computed 779/779; see tex).
Enumerate: T_L in {B,C}, T_E = other; E_0 = any good T_E-configuration other than (1,a_{T_E}); L_0 = (1,a_{T_L}),
translated so that the squares coincide; t_E1, t_L1 = the type of the next step at the created rhombus.
Constraints (tags) -- see the tex, Lemma loc2:
  M0   : every tile known to be in M0 -- H_E, H_L, G(E_1) minus C(E_0), G(L_1) minus C(L_0) -- pairwise equal or
         interior-disjoint (M0 is a tiling).
  BRDE : G(E_1) lies on X = M0 - H_E + C(E_0): each tile is in C(E_0), or it is an M0 tile not in H_E; and the tiles of
         X that we know (C(E_0), M0-known minus H_E) are pairwise compatible.  BRDL symmetric.
  UQE  : (Lemma U at R_E^1) no type t != t_E1 has G_t(R_E^1) inside the known tiles of X.  UQL symmetric.
  UQ0  : (Lemma U at the contacts) no type t != t_E0 has G_t(R_E) inside M0-known; same for L.
LESION=<tag> drops one constraint.  Output: survivors, and (Share_2) on each."""
import os, sys, itertools
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
LES = set(os.environ.get('LESION', '').split(',')) - {''}
aT = {'B': 'a1', 'C': 'a0'}
def sqof(S): return next(t for t in S if t[0] == 'sq')
def pairwise_ok(S): return all(a == b or not overlap(a, b) for a, b in itertools.combinations(list(S), 2))
def configs_at(T, R):
    out = {}
    for typ in TYPES:
        if (T, classkey(R), typ) in CF: out[typ] = step(T, R, typ)
    return out
surv = []; st = Counter()
for TL in 'BC':
    TE = OTHER[TL]
    for (T_, ckE, tE0), (RcE, remE, RnE, creE) in CF.items():
        if T_ != TE: continue
        if cls(TE, RcE) == 1 and tE0 == aT[TE]: continue          # Prop orient: E is not at (1, a_TE)
        RE = RcE; HE = set(remE) | {RE}; CE0 = set(creE) | {RnE}
        u = sqof(remE)
        RcL, remL, RnL, creL = CF[(TL, K1[TL], aT[TL])]
        uL = sqof(remL)
        if shape_ok := (sorted(tuple(a - b for a, b in zip(P, base(u))) for P in u[1]) ==
                        sorted(tuple(a - b for a, b in zip(P, base(uL))) for P in uL[1])):
            pass
        if not shape_ok: continue
        D = tuple(a - b for a, b in zip(base(u), base(uL)))
        RL = tr(RcL, D); HL = {tr(t, D) for t in remL} | {RL}; CL0 = {tr(t, D) for t in creL} | {tr(RnL, D)}
        cE0 = (cls(TE, RE), tE0); st[('contact pairs (shape match)',)] += 1
        if not pairwise_ok(set(nonrh(HE)) | set(nonrh(HL))): st[('contact pair killed by M0',)] += 1; continue
        RE1 = tr(RnE, (0, 0, 0, 0)); RL1 = tr(RnL, D)
        cfE, cfL = configs_at(TE, RE1), configs_at(TL, RL1)
        for tE1, tL1 in itertools.product(cfE, cfL):
            remE1, _, creE1 = cfE[tE1]; remL1, _, creL1 = cfL[tL1]
            GE1, GL1 = set(nonrh(remE1)), set(nonrh(remL1))
            nHE, nHL = set(nonrh(HE)), set(nonrh(HL))
            fail = []
            # NR: no-revisit (REM+CRE+MU) -- G(E_1) is not created by E_0, and an M0 tile is removed once
            if (GE1 & CE0) or (GE1 & nHE) or (GL1 & CL0) or (GL1 & nHL): fail.append('NR')
            # M0: all non-rhombus tiles removed by J_E or J_L are M0 tiles
            if not pairwise_ok(nHE | nHL | GE1 | GL1): fail.append('M0')
            # XE0: board of J_E before p_E holds H_E, and the W_L tiles nHL, GL1 (minimality of p_E), and GE1
            if not pairwise_ok(HE | nHL | GL1 | GE1): fail.append('XE0')
            if not pairwise_ok(HL | nHE | GE1 | GL1): fail.append('XL0')
            XE = CE0 | ((nHL | GL1 | GE1) - HE); XL = CL0 | ((nHE | GE1 | GL1) - HL)
            if not pairwise_ok(XE): fail.append('XE1')
            if not pairwise_ok(XL): fail.append('XL1')
            if any(t != tE1 and set(nonrh(cfE[t][0])) <= XE for t in cfE): fail.append('UQE')
            if any(t != tL1 and set(nonrh(cfL[t][0])) <= XL for t in cfL): fail.append('UQL')
            c0E, c0L = configs_at(TE, RE), configs_at(TL, RL)
            if any(t != tE0 and set(nonrh(c0E[t][0])) <= (HE | nHL | GL1 | GE1) for t in c0E) or \
               any(t != aT[TL] and set(nonrh(c0L[t][0])) <= (HL | nHE | GE1 | GL1) for t in c0L): fail.append('UQ0')
            live = [f for f in fail if f not in LES]
            st[('combo', 'killed by', tuple(fail))] += 1
            if live: continue
            OE, OL = GE1 - nHE - nHL, GL1 - nHE - nHL
            s2 = (OE == OL, len(OE), len(GE1 & HL), len(GL1 & HE), bool(GE1 & CE0), bool(GL1 & CL0),
                  len((set(creE1) | {cfE[tE1][1]}) & (set(creL1) | {cfL[tL1][1]})))
            surv.append(dict(TL=TL, cE0=cE0, tE1=tE1, tL1=tL1, s2=s2, fail=fail))
print('LESION', sorted(LES) or 'none')
for k in sorted(st, key=str): print(st[k], k)
print('survivors', len(surv))
for s in surv:
    print('SURV T_L', s['TL'], 'E0', s['cE0'], 'E1', s['tE1'], 'L1', s['tL1'],
          '| Oth equal, |Oth|, |G(E1)&H_L|, |G(L1)&H_E|, G(E1) meets C(E0), G(L1) meets C(L0), |C(E1)&C(L1)| =', s['s2'],
          '| lesioned fails', s['fail'])
print('SHARE2 on all survivors:', all(s['s2'][0] and s['s2'][1] == 2 for s in surv))
