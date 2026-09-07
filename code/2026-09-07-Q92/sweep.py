import sys, itertools, sympy as sp
sys.path.insert(0, '/home/clio/projects/scratch/q92')
import engineA as A, engineB as B
t = sp.Symbol('t')

ok = bad = 0; failures = []
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
EMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5
pairs = [(e, f) for e in range(1, EMAX+1) for f in range(1, EMAX+1) if e < f]
lams = [l for n in range(0, NMAX+1) for l in A.partitions(n)]
for (e, f) in pairs:
    for lam in lams:
        a = A.commutator(lam, e, f)
        b = B.closed_form(lam, e, f)
        keys = set(a) | set(b)
        agree = all(sp.simplify(sp.expand(a.get(k, 0) - b.get(k, 0))) == 0 for k in keys)
        if agree: ok += 1
        else:
            bad += 1; failures.append((e, f, lam, a, b))
print(f"AGREE {ok} / {ok+bad}   (|lam|<={NMAX}, 1<=e<f<={EMAX})")
for x in failures[:6]: print("FAIL", x)
