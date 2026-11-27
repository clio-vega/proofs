"""R3: tiles removed by E^L_{k+2}, E^L_{k+3} that are not in C(L_{k+2}), G(E_{k+2}), G(E_{k+3}), C(E^L_{k+2}).
For each: is it created by some L_j (corridor, j>=2) or L_0,L_1 (then not M0: harmless)?  else is g+alpha in G(E_j)
for which j (AW-shadow)?  Over all words up to length 6 and all contacts."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
H = Counter()
for TL in 'BC':
    TE = OTHER[TL]; al = ALPHA[TL]; mal = tuple(-a for a in al)
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, TE):
        if 'b' not in configs_at(TL, RL1): continue
        remL1, R2, creL1 = configs_at(TL, RL1)['b']
        for k in range(0, 7):
            for word in itertools.product(('a0', 'a1'), repeat=k):
                R = R2; GE = {}; CLs = {'L0': set(nonrh(CL0)), 'L1': set(nonrh(creL1))}
                for j, t in enumerate(word, start=2):
                    rem, Rn, cre = step(TE, R, t); GE[j] = (R, set(nonrh(rem))); CLs['L%d' % j] = set(nonrh(cre)); R = Rn
                Rk2 = R
                remEb, RE3, creEb = step(TE, Rk2, 'b'); GE['k+2'] = (Rk2, set(nonrh(remEb)))
                remLb, _, creLb = step(TL, Rk2, 'b'); CLs['Lk+2'] = set(nonrh(creLb))
                rem2, R3, cre2 = step(TE, tr(Rk2, mal), aT[TL])
                for x, y in (('b', aT[TE]), (aT[TL], 'b')):
                    rx = step(TE, RE3, x)[0]; GE['k+3'] = (RE3, set(nonrh(rx)))
                    ry = step(TE, R3, y)[0]
                    lft = (set(nonrh(rem2)) | set(nonrh(ry))) - CLs['Lk+2'] - GE['k+2'][1] - GE['k+3'][1] - set(cre2)
                    for g in lft:
                        cr = [n for n, S in CLs.items() if g in S]
                        sh = [j for j, (_, S) in GE.items() if tr(g, al) in S]
                        H[(TL, x, g[0], 'created by', tuple(cr), 'g+alpha in G(E_j)', tuple(sh) if k == 0 else tuple(str(s).replace(str(k+1), 'k+1') for s in sh))] += 1
for kx in sorted(H, key=str): print(H[kx], kx)
