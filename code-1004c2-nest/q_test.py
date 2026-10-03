import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from onesided import F
from gen import is_pf2
from itertools import product
from collections import Counter

print('=== (Q) exhaustive, small ===', flush=True)
for m in (2, 3, 4):
    st = Counter(); wit = []
    Um = {2: 7, 3: 5, 4: 4}[m]
    for L in product(range(1, 3), repeat=m):
        for U in product(*[range(L[i], Um+1) for i in range(m)]):
            for T in range(sum(L), sum(U)+1):
                f, nW = F(list(L), list(U), T)
                if f is None: continue
                st['cases'] += 1
                st['maxW'] = max(st.get('maxW', 0), nW)
                st['maxdeg'] = max(st.get('maxdeg', 0), len(f)-1)
                if is_pf2(f): st['PF2'] += 1
                else:
                    st['FAIL'] += 1
                    if len(wit) < 4: wit.append((L, U, T, nW, f))
    print(' m=%d' % m, dict(st), flush=True)
    for w in wit: print('    FAIL', w, flush=True)

print('=== (Q) random, wide ===', flush=True)
rng = random.Random(7); st = Counter(); wit = []
for _ in range(60000):
    m = rng.randint(2, 6)
    L = [rng.randint(1, 4) for _ in range(m)]
    U = [L[i] + rng.randint(0, 9) for i in range(m)]
    T = rng.randint(sum(L), sum(U))
    f, nW = F(L, U, T)
    if f is None or nW == 0: continue
    st['cases'] += 1
    st['maxW'] = max(st.get('maxW', 0), nW)
    st['maxdeg'] = max(st.get('maxdeg', 0), len(f)-1)
    if is_pf2(f): st['PF2'] += 1
    else:
        st['FAIL'] += 1
        if len(wit) < 4: wit.append((m, L, U, T, nW, f))
print(dict(st), flush=True)
for w in wit: print('    FAIL', w, flush=True)
