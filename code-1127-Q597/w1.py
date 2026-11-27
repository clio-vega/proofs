"""W1 (for Q2'): no non-rhombus tile removed by E_{k+3} (any of the 3 types at class 1) is removed by
L_0, L_1 (contact), by E^L_j = E_j - alpha (2<=j<=k+1), or by E^L_{k+2} = T_E's (0,a_{T_L}) at R_{k+2}-alpha.
Over every case-b contact class and every corridor word t_2..t_{k+1} in {a0,a1}^k; window phi bounds k.
Also W1x: same for the tiles removed by E^L_{k+3} (both candidates) vs E_j, for the symmetric use (instrument)."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
PHI = (-1, -2, 1, -2)
f = lambda P: sum(p * x for p, x in zip(PHI, P))
mn = lambda R: min(R[1])
def rng(tiles, R):
    v = [f(P) - f(mn(R)) for g in tiles for P in g[1]]; return min(v), max(v)
def run(TL, word, RL1, HL, track):
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL])
    remL1, R2, creL1 = configs_at(TL, RL1)['b']
    forb = {}  # tile -> who removes it (on the L / E^L side)
    for g in nonrh(HL): forb.setdefault(g, 'L0')
    for g in nonrh(remL1): forb.setdefault(g, 'L1')
    R = R2
    for j, t in enumerate(word, start=2):
        assert cls(TE, R) == 0
        remEj, Rn, _ = step(TE, R, t)
        remELj = step(TE, tr(R, mal), t)[0]
        for g in nonrh(remELj): forb.setdefault(g, 'EL%d' % j)
        R = Rn
    Rk2 = R; assert cls(TE, Rk2) == 0
    for g in nonrh(step(TE, tr(Rk2, mal), aT[TL])[0]): forb.setdefault(g, 'ELk+2')
    remEb, RE3, _ = step(TE, Rk2, 'b')
    hits = []
    for x, (rx, _, _) in configs_at(TE, RE3).items():
        for h in nonrh(rx):
            if h in forb: hits.append((x, h[0], forb[h]))
        track.append(rng(nonrh(rx), Rk2))
    return hits, R2
# windows
trE, trL = [], []
for TL in 'BC':
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, OTHER[TL]):
        if 'b' not in configs_at(TL, RL1): continue
        _, R2 = run(TL, (), RL1, HL, trE)
        remL1 = configs_at(TL, RL1)['b'][0]
        trL.append(rng(nonrh(HL) + nonrh(remL1), R2))
a1 = min(r[0] for r in trE); a2 = max(r[1] for r in trE); b1 = min(r[0] for r in trL); b2 = max(r[1] for r in trL)
amax = max(abs(f(a)) for a in ALPHA.values())
# E^L_j tiles: vertices of removed tiles of a class-0 a-step at R have phi - phi(min R) in [-4,1] (persist.py), shifted by -alpha
NL = b2 - a1          # need k + a1 <= b2
NEL = 1 + amax - a1   # need (k+2-j) + a1 <= 1 + amax
print('E_{k+3} tiles phi-range rel R_{k+2}: [%d,%d]; L0,L1 tiles rel R_2: [%d,%d]; bound k<=%d (L), k+2-j<=%d (E^L)' % (a1, a2, b1, b2, NL, NEL))
N = max(NL, NEL) + 1
tot = 0; H = Counter()
for TL in 'BC':
    for cE0, RE, HE, CE0, RE1, RL, HL, CL0, RL1, u in contacts(TL, OTHER[TL]):
        if 'b' not in configs_at(TL, RL1): continue
        for k in range(0, N + 1):
            for word in itertools.product(('a0', 'a1'), repeat=k):
                hits, _ = run(TL, word, RL1, HL, [])
                tot += 1
                for hh in hits: H[(TL, k) + hh] += 1
print('W1 words checked', tot, '(k up to %d)' % N, 'hits:', dict(H) if H else 'NONE')
