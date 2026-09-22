from itertools import combinations
from beads import *
import collections
def all_beadsets(n):
    return [frozenset(c) for k in range(1, n) for c in combinations(range(n), k)]
def hops(A,n): return sum(1 for p in A if (p+1)%n not in A)

# non-degenerate witnesses
w1 = w2 = None
for n in range(3, 8):
    for S, T, win in additive_pairs(n):
        if not S or not T: continue
        E = clio_letters(S, T, n); hs = {}
        for A in all_beadsets(n):
            nu = act_word(A, word_cd(T,n), n)
            if nu is None: continue
            if act_word(nu, word_cd(S,n), n) is None: continue
            hs[A] = nu
        if not hs: continue
        vals = {hops(nu,n) for nu in hs.values()}
        if len(vals) > 1 and w1 is None:
            w1 = (n, sorted(S), sorted(T), [(sorted(A), sorted(nu), hops(nu,n)) for A,nu in hs.items()], len(E))
        if w2 is None:
            for A,nu in hs.items():
                if len(E) != hops(nu,n):
                    w2 = (n, sorted(S), sorted(T), sorted(A), sorted(nu), len(E), hops(nu,n)); break
    if w1 and w2: break
print("K1 smallest non-degenerate witness (#hops(nu) not an invariant of v):")
print("   n=%d S=%s T=%s  #moves=%d" % (w1[0], w1[1], w1[2], w1[4]))
for A,nu,h in w1[3]: print("      A=%s -> nu=%s  hops=%d" % (A,nu,h))
print("K2 smallest non-degenerate (C)-failure: n=%d S=%s T=%s A=%s nu=%s  #moves=%d vs hops(nu)=%d" % w2)

# Is the arc multiset S + T an invariant of v across the WHOLE fibre?
print()
print("Is the multiset S + T constant on the fibre F(v)?")
c = collections.Counter()
for n in range(3, 8):
    fib = collections.defaultdict(list)
    for S, T, win in additive_pairs(n):
        fib[win].append((S, T))
    for win, lst in fib.items():
        ms = {tuple(sorted(collections.Counter(list(S)+list(T)).items())) for S, T in lst}
        c['fibres'] += 1
        c['constant'] += (len(ms) == 1)
    print(f"   n<={n}: fibres {c['fibres']}, arc-multiset constant on {c['constant']}")

print()
print("Restricted to 321-avoiding v (= realisable fibres):")
def is321(win, n):
    def W(i):
        q, r = divmod(i-1, n); return win[r] + q*n
    hi = 3*n
    for x in range(1, n+1):
        for y in range(x+1, hi):
            if W(y) >= W(x): continue
            for z in range(y+1, hi):
                if W(z) < W(y): return False
    return True
c2 = collections.Counter(); ex = None
for n in range(3, 8):
    fib = collections.defaultdict(list)
    for S, T, win in additive_pairs(n):
        fib[win].append((S, T))
    for win, lst in fib.items():
        if not is321(win, n): continue
        c2['fib321'] += 1
        ms = {tuple(sorted(collections.Counter(list(S)+list(T)).items())) for S, T in lst}
        c2['const'] += (len(ms) == 1)
        if len(ms) > 1 and ex is None:
            ex = (n, win, [(sorted(S), sorted(T)) for S, T in lst])
    print(f"   n<={n}: 321-avoiding fibres {c2['fib321']}, arc-multiset constant on {c2['const']}")
print("   first non-constant 321 fibre:", ex)
