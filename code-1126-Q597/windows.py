"""Distinct front/back windows up to one common Z^4 translation (origin = min point of the shared contact square)."""
import pickle, glob, sys
from collections import Counter
sys.argv = ['x', '5', '1', '24']
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
def shift(s, o):
    f = lambda S: frozenset((t[0], frozenset(tuple(a-b for a, b in zip(p, o)) for p in t[1])) for t in S)
    return (s[0], f(s[1]), f(s[2]))
st = Counter(); FR = Counter(); BKc = Counter()
for fn in sorted(glob.glob('/tmp/step0_5_*.pkl')):
    for case, L, info in pickle.load(open(fn, 'rb')):
        E = 1 - L
        TE, AE, ArE, pE, laE, lbE = info[E]; TL, AL, ArL, pL, laL, lbL = info[L]
        sq = [t for t in AL[pL][1] if t[0] == 'sq']
        u = [t for t in AE[pE][1] if t in sq]
        o = min(u[0][1])
        if case == 'b':
            k = laL - 3
            front = tuple(shift(s, o) for s in (AE[pE], AE[pE+1], ArE[pE], ArE[pE+1], AL[pL], AL[pL+1], ArL[pL]))
            FR[(TL, front)] += 1
            ob = min(p for t in AE[pE+k+2][1] for p in t[1])   # back origin: E's first post-corridor step
            back = tuple(shift(s, ob) for s in (AE[pE+k+2], AE[pE+k+3], ArE[pE+k+2], ArE[pE+k+3], AL[pL+k+2], ArL[pL+k+1], ArL[pL+k+2]))
            BKc[(TL, back)] += 1
        else:
            front = tuple(shift(s, o) for s in list(AE[pE:pE+3]) + list(ArE[pE:pE+3]) + list(AL[pL:pL+2]) + list(ArL[pL:pL+2]))
            FR[('a', TL, front)] += 1
print('case b: distinct front windows', len([k for k in FR if k[0] != 'a']), 'over', sum(v for k, v in FR.items() if k[0] != 'a'))
print('case b: distinct back windows', len(BKc), 'over', sum(BKc.values()))
print('case a: distinct whole-middle windows', len([k for k in FR if k[0] == 'a']), 'over', sum(v for k, v in FR.items() if k[0] == 'a'))
print('front multiplicities', sorted([v for k, v in FR.items() if k[0] != 'a'], reverse=True))
print('back multiplicities', sorted(BKc.values(), reverse=True))
print('case a multiplicities', sorted([v for k, v in FR.items() if k[0] == 'a'], reverse=True))
