"""Pair lemma on E^L's board before E^L_{k+3}.  Present there (proved/claimed):
  P1 = G(E_{k+3}) minus G(E^L_{k+2})   [Q2'': W1 window shows only E^L_{k+2} removes a tile of G(E_{k+3})]
  P2 = nonrh C(E^L_{k+2})               [just created]
  P3 = G(E^L_{k+3})                     [Lemma word]
For each x (type of E_{k+3}, class 1) and y in {a_{T_E}, b} (Prop wb(1)): pairwise compatible?  Census pairs must be the only
compatible ones.  Lesion: drop P1 (only P2+P3): should NOT determine the pair."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
def compat(S): return not killers(set(S), set(S)) if False else all(a == b or not overlap(a, b) for a, b in itertools.combinations(list(S), 2))
ok = True
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remEb, RE3, creEb = step(TE, R, 'b')
    rem2, R3, cre2 = step(TE, Ra, aT[TL]); G2 = set(nonrh(rem2))
    cE = configs_at(TE, RE3); cEL = configs_at(TE, R3)
    for x, (rx, _, _) in cE.items():
        for y in (aT[TE], 'b'):
            P1 = set(nonrh(rx)) - G2; P2 = set(nonrh(cre2)); P3 = set(nonrh(cEL[y][0]))
            c = compat(P1 | P2 | P3); cl = compat(P2 | P3)
            good = (x, y) in {('b', aT[TE]), (aT[TL], 'b')}
            print('T_L', TL, 'x', x, 'y', y, 'census' if good else 'bad', '| compatible:', c, '| lesion (no P1):', cl,
                  '| |G_x & G(EL_k+2)| =', len(set(nonrh(rx)) & G2))
            ok &= (c == good)
print('PAIR LEMMA on E^L board before E^L_{k+3}:', ok)
print('--- lesion 2: the uncorrected premise (ALL of G(E_{k+3}) present) together with C(E^L_{k+2})')
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remEb, RE3, creEb = step(TE, R, 'b'); rem2, R3, cre2 = step(TE, Ra, aT[TL])
    cE = configs_at(TE, RE3); cEL = configs_at(TE, R3)
    print('T_L', TL, 'compatible pairs:', [(x, y) for x in cE for y in (aT[TE], 'b')
          if compat(set(nonrh(cE[x][0])) | set(nonrh(cre2)) | set(nonrh(cEL[y][0])))])
