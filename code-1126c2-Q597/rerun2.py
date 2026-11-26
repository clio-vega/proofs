"""PROVE 11-26 c2: E's re-run front E^L_0, E^L_1 by presence + Lemma U, on the 8 local classes of local2.py.
Known-present tiles on E^L's board Y_p (before step p_E), all proved facts (see tex, Prop rerun):
   K = nonrh(C(L_0) u C(L_1))            -- created by J_L, final (Lemma fin), not M0 (no re-creation) => in Xadd_L c Y_p (dinv)
     u (nonrh(H_E) - nonrh(H_L) - G(L_1)) -- M0 tiles on X_p not removed by J_L (prefix: minimality; corridor/exit = E's
                                              steps, which remove them only at E_0; suffix = L^E's steps on M0.J_E, where
                                              H_E tiles are absent)
   plus R_E (the re-run's rhombus at p is E's, common prefix).
Lemma U: exactly one good type is present at R_E on Y_p, and it is the step taken (E^L steps at p: Thm fc).
So if exactly one type t' has G_{t'}(R_E) c K, then E^L_0 = t'.  Then iterate once at E^L_0's new rhombus with
   Kn = (K - G(E^L_0)) u nonrh(C(E^L_0)).
Also prints compatibility of K (must be pairwise compatible: subset of the tiling Y_p) and compares with the census."""
import itertools, pickle
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
exec(open('/home/clio/projects/proofs/code-1126c2-Q597/why2.py').read().split("MODE = sys.argv")[0].split('"""', 2)[2])
def pairwise_ok(S): return all(a == b or not overlap(a, b) for a, b in itertools.combinations(list(S), 2))
def configs_at(T, R): return {typ: step(T, R, typ) for typ in TYPES if (T, classkey(R), typ) in CF}
OUT = []
for TL in 'BC':
    TE = OTHER[TL]
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        nHE, nHL = set(nonrh(HE)), set(nonrh(HL))
        cfE, cfL = configs_at(TE, RE1), configs_at(TL, RL1)
        for tE1, tL1 in itertools.product(cfE, cfL):
            GE1, GL1 = set(nonrh(cfE[tE1][0])), set(nonrh(cfL[tL1][0]))
            if not pairwise_ok(nHE | nHL | GE1 | GL1): continue
            caseb = (tL1 == 'b')
            CL1 = set(cfL[tL1][2])
            K = set(nonrh(CL0 | CL1)) | (nHE - nHL - GL1)
            okK = pairwise_ok(K | {RE})
            c0 = configs_at(TE, RE); tE0 = cE0[1]
            pres = [t for t in c0 if t != tE0 and set(nonrh(c0[t][0])) <= K]
            line = dict(TL=TL, E0=cE0, E1=tE1, L1=tL1, case='b' if caseb else 'a', Kcompat=okK, present0=pres)
            if len(pres) == 1:
                t0 = pres[0]; rem0, R1, cre0 = c0[t0]
                Kn = (K - set(nonrh(rem0))) | set(nonrh(cre0))
                c1 = configs_at(TE, R1)
                pres1 = [t for t in c1 if set(nonrh(c1[t][0])) <= Kn]
                line.update(EL0=(cls(TE, RE), t0), EL1_class=cls(TE, R1), present1=pres1, K1compat=pairwise_ok(Kn | {R1}))
                # source of E^L_0's tiles
                line['EL0 tiles from'] = sorted(('C(L0)' if g in CL0 else 'C(L1)' if g in CL1 else 'H_E' if g in nHE else '?') for g in nonrh(rem0))
                if len(pres1) == 1:
                    rem1 = c1[pres1[0]][0]
                    line['EL1 tiles from'] = sorted(('C(L0)' if g in CL0 else 'C(L1)' if g in CL1 else 'C(EL0)' if g in cre0 else 'H_E' if g in nHE else '?') for g in nonrh(rem1))
            OUT.append(line)
for l in OUT: print(l)
# ---- census comparison (exact tile sets in the shared-square frame); rerun_census.pkl from the inline census script
import os
LOCAL = set()
for TL in 'BC':
    TE = OTHER[TL]
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        nHE, nHL = set(nonrh(HE)), set(nonrh(HL))
        for tL1 in ('b', {'B': 'a0', 'C': 'a1'}[TL] if False else [t for t in configs_at(TL, RL1) if t not in ('b', {'B':'a1','C':'a0'}[TL])][0]):
            rL1 = configs_at(TL, RL1)[tL1]; GL1 = set(nonrh(rL1[0])); CL1 = set(rL1[2])
            K = set(nonrh(CL0 | CL1)) | (nHE - nHL - GL1)
            c0 = configs_at(TE, RE); pres = [t for t in c0 if set(nonrh(c0[t][0])) <= K]
            assert len(pres) == 1, pres
            rem0, R1, cre0 = c0[pres[0]]; Kn = (K - set(nonrh(rem0))) | set(nonrh(cre0))
            c1 = configs_at(TE, R1); pres1 = [t for t in c1 if set(nonrh(c1[t][0])) <= Kn]; assert len(pres1) == 1
            b = min(u[1]); f = lambda S: frozenset((t[0], frozenset(tuple(a - c for a, c in zip(P, b)) for P in t[1])) for t in S)
            LOCAL.add((TL, f(HE), f(set(rem0)), f(set(c1[pres1[0]][0])), f(set(rL1[0]) )))
fn = '/home/clio/projects/proofs/code-1126c2-Q597/rerun_census.pkl'
if os.path.exists(fn):
    C = pickle.load(open(fn, 'rb'))
    print('CENSUS re-run front classes', len(C), 'local', len(LOCAL), 'census in local', sum(k in LOCAL for k in C), 'local seen', len(LOCAL & set(C)))
