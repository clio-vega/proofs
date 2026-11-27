"""After the back window: do E_{k+3} (on M0) and E^L_{k+3} create the same rhombus, for the two census pairs?
Also the local part of D = X* \\ Y* that J_E could read: M0 tiles removed by E^L's back steps but not by E_{k+2},E_{k+3}."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remEb, RE3, creEb = step(TE, R, 'b')
    rem2, R3, cre2 = step(TE, Ra, aT[TL])
    for x, y in (('b', aT[TE]), (aT[TL], 'b')):
        rx, RE4, cx = step(TE, RE3, x); ry, REL4, cy = step(TE, R3, y)
        print('T_L', TL, 'pair', x, y, '| same rhombus after:', RE4 == REL4, '| class', cls(TE, RE4),
              '| created sets equal:', set(cx) == set(cy), '| G(E^L_{k+2..k+3}) - G(E_{k+2..k+3}) - C(E^L_{k+2}):',
              sorted(g[0] for g in (set(nonrh(rem2)) | set(nonrh(ry))) - set(nonrh(remEb)) - set(nonrh(rx)) - set(cre2)))
