"""PROVE 11-21 c2: base-point displacement (Z^4) of every good configuration, per target, class, type; and a search for a
linear functional phi on Z^4 with phi >= 0 on all T-steps and phi <= 0 on all T'-steps (a separator of the two futures
after the split).  Class-0 a-steps are shared (Lemma shared), so phi must vanish on them.   usage: python3 stepdisp.py"""
import os, sys, itertools
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '.')
exec(open('words.py').read().split('\ndef run(T, L, st):')[0])
D = {}
for T in 'BC':
    for (ck, typ), (Rc, rem, Rn, cre) in configs(T).items():
        d = tuple(b - a for a, b in zip(base(Rc), base(Rn)))
        D[(T, tuple(sorted(ck))[0], typ)] = d
        print(T, sorted(ck), typ, 'base displacement', tuple(str(x) for x in d))
best = []
R = range(-3, 4)
for phi in itertools.product(R, repeat=4):
    if not any(phi): continue
    val = lambda d: sum(Fr(p) * x for p, x in zip(phi, d))
    okB = all(val(d) >= 0 for (T, c, t), d in D.items() if T == 'B')
    okC = all(val(d) <= 0 for (T, c, t), d in D.items() if T == 'C')
    if okB and okC:
        strictB = sum(1 for (T, c, t), d in D.items() if T == 'B' and val(d) > 0)
        strictC = sum(1 for (T, c, t), d in D.items() if T == 'C' and val(d) < 0)
        best.append((strictB + strictC, phi, strictB, strictC))
best.sort(reverse=True)
print('separators found', len(best)); print(best[:8])
