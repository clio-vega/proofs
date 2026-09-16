from twobead import words
from supports import enumerate_supports
from itertools import product

def Cset(d, delta, b):
    return frozenset(w for w in words(d-1) if w[delta-1]==b and w[d-delta-1]==1-b)

def all_C(d):
    out = {}
    for delta in range(1, d):
        if 2*delta == d: continue
        for b in (0,1):
            out[(delta,b)] = Cset(d,delta,b)
    return out

for d in range(3, 8):
    ws, good = enumerate_supports(d)
    prop = [Z for Z in good if 0 < len(Z) < len(ws)]
    Cs = all_C(d)
    distinct = set(Cs.values())
    contained = [Z for Z in prop if any(Z <= C for C in distinct)]
    print('d=%d: %d proper; %d distinct C-sets (size %d); contained in some C: %d/%d'
          % (d, len(prop), len(distinct), 2**(d-3), len(contained), len(prop)))
    if len(contained) != len(prop):
        bad = [Z for Z in prop if not any(Z<=C for C in distinct)]
        print('    NOT CONTAINED:', [sorted(''.join(map(str,w)) for w in Z) for Z in bad[:5]])
    # singleton criterion
    singles = sorted(w for Z in prop if len(Z)==1 for w in Z)
    n1 = [w for w in words(d-1) if w[0]!=w[d-2]]
    print('    singletons admissible: %d ; words with w_1 != w_{d-1}: %d ; subset? %s'
          % (len(singles), len(n1), set(singles) <= set(n1)))
    missing = [w for w in n1 if w not in set(singles)]
    if missing: print('    w_1!=w_{d-1} but NOT admissible:', [''.join(map(str,w)) for w in missing])
