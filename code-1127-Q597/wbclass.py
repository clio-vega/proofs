"""(W_B) from the local derivation: build the back window (E_{k+2},E_{k+3},E^L_{k+2},E^L_{k+3},L_{k+2},L^E_{k+1},L^E_{k+2})
from the proved local facts, for each T_L and each of the two pairs, and compare with the 4 census representatives
in code-1126-Q597/windows.pkl (translated so E_{k+2}'s least point is 0).  Then (Bk) holds by Lemma win."""
import pickle
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
aT = {'B': 'a1', 'C': 'a0'}
W = pickle.load(open('/home/clio/projects/proofs/code-1126-Q597/windows.pkl', 'rb'))
reps = {key: v for key, v in W['back'].items()}
s0 = next(iter(reps))[1][0]
print('census step format:', s0[0], type(s0[1]), sorted(t[0] for t in s0[1]), sorted(t[0] for t in s0[2]))
def conv(typ, R, res):
    rem, Rn, cre = res
    F = lambda S: frozenset((t[0], frozenset(t[1])) for t in S)
    g = set(rem) | {R}; n = set(cre) | {Rn}
    return (typ, F(g), F(n))
def shift(s, o):
    f = lambda S: frozenset((t[0], frozenset(tuple(a - b for a, b in zip(p, o)) for p in t[1])) for t in S)
    return (s[0], f(s[1]), f(s[2]))
found = 0
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    EB = step(TE, R, 'b'); RE3 = EB[1]
    LB = step(TL, R, 'b'); RL3 = LB[1]
    EL2 = step(TE, Ra, aT[TL]); R3 = EL2[1]
    LE1 = step(TL, Ra, 'b'); LE2 = step(TL, tr(RL3, mal), aT[TL])
    for x, y in (('b', aT[TE]), (aT[TL], 'b')):
        win = [conv('b', R, EB), conv(x[0], RE3, step(TE, RE3, x)), conv('a', Ra, EL2), conv(y[0], R3, step(TE, R3, y)),
               conv('b', R, LB), conv('b', Ra, LE1), conv('a', tr(RL3, mal), LE2)]
        o = min(p for t in win[0][1] for p in t[1])
        win = tuple(shift(s, o) for s in win)
        hit = [(key[0], v[1]) for key, v in reps.items() if key[0] == TL and key[1] == win]
        hit2 = [(key[0], v[1]) for key, v in reps.items() if key[0] == TL and tuple(s[1:] for s in key[1]) == tuple(s[1:] for s in win)]
        print('T_L', TL, 'pair', x, y, '| matches census rep (exact):', hit, '| ignoring typ field:', hit2)
        found += bool(hit or hit2)
print('local windows matched:', found, 'of 4; census back reps:', len(reps))
