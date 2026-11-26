"""Base of Merge': is E^L_1 a translate of L_1 (case b)?  If so, alpha := that translation, and R_2 - alpha is created by E^L_1."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c2-Q597/rerun2.py').read().split("OUT = []")[0])
def shape(S):
    b = min(P for t in S for P in t[1]); return frozenset((t[0], tuple(sorted(tuple(a-c for a, c in zip(P, b)) for P in t[1]))) for t in S), b
ALPHA = {}
for TL in 'BC':
    TE = OTHER[TL]
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        cfL = configs_at(TL, RL1)
        if 'b' not in cfL: continue
        remL1, RL2, creL1 = cfL['b']
        nHE, nHL = set(nonrh(HE)), set(nonrh(HL)); GL1 = set(nonrh(remL1)); CL1 = set(creL1)
        K = set(nonrh(CL0 | CL1)) | (nHE - nHL - GL1)
        c0 = configs_at(TE, RE); pres = [t for t in c0 if set(nonrh(c0[t][0])) <= K]; assert len(pres) == 1
        rem0, R1, cre0 = c0[pres[0]]; Kn = (K - set(nonrh(rem0))) | set(nonrh(cre0))
        c1 = configs_at(TE, R1); pres1 = [t for t in c1 if set(nonrh(c1[t][0])) <= Kn]; assert pres1 == ['b']
        remE1, RE2L, creE1 = c1['b']
        sL, bL = shape(set(remL1) | {RL1} | set(creL1) | {RL2}); sE, bE = shape(set(remE1) | {R1} | set(creE1) | {RE2L})
        alpha = tuple(a - b for a, b in zip(bL, bE))
        print(TL, 'E0', cE0, 'E^L_1 translate of L_1 (removed+created+rhombi):', sL == sE, 'alpha (L - E^L)', alpha,
              '| created rhombus of E^L_1 == R2 - alpha:', tr(RL2, tuple(-a for a in alpha)) == RE2L)
        ALPHA.setdefault(TL, set()).add(alpha)
print('ALPHA by T_L', ALPHA)
