"""(i) packed_at shape: the rhombus {q,q+u_T,q+w_T,q+u_T+w_T} has classkey K1[T], not K0 (so no journey stops at class 0).
(ii) persistence: for corridor words t_j..t_i (class-0 a-steps, each creating the next rhombus), with i>j>=2,
     G_{t_j}(R_j - alpha) (E^L_j, by induction) never contains a tile of nonrh C(L_i) = nonrh C_{t_i}(R_i).
     Monotone window: phi linear with phi(d_a0), phi(d_a1) > 0 bounds i-j; all words up to that length are checked."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
U = {'B': (-1, 1, 0, 0), 'C': (1, 0, 0, 0)}; W = {'B': (0, 0, 1, -1), 'C': (0, 0, -1, 0)}
for T in 'BC':
    q = (0, 0, 0, 0); u, w = U[T], W[T]
    R = ('rh', tuple(sorted([q, u, w, tuple(a + b for a, b in zip(u, w))])))
    print('(i)', T, 'packed shape classkey == K1[T]:', classkey(R) == K1[T], '== K0:', classkey(R) == K0, '| other target class1:', classkey(R) == K1[OTHER[T]])
add = lambda P, v: tuple(a + b for a, b in zip(P, v))
def mn(R): return min(R[1])
R0 = canon('B', 0)
D = {}
for t in ('a0', 'a1'):
    rem, Rn, cre = step('B', R0, t); D[t] = tuple(a - b for a, b in zip(mn(Rn), mn(R0)))
    assert cls('B', Rn) == 0
print('displacements of min R at class 0:', D)
# search a small integer phi positive on both
best = None
for phi in itertools.product(range(-2, 3), repeat=4):
    v = [sum(p * d for p, d in zip(phi, D[t])) for t in D]
    if min(v) > 0: best = phi; break
print('phi', best, [sum(p * d for p, d in zip(best, D[t])) for t in D])
f = lambda P: sum(p * x for p, x in zip(best, P))
def walk(word, TE, R):
    out = []
    for t in word:
        rem, Rn, cre = step(TE, R, t); out.append((R, set(nonrh(rem)), set(nonrh(cre)))); R = Rn
    return out
# windows
lo = hi = None
for T in 'BC':
    for t in ('a0', 'a1'):
        rem, Rn, cre = step(T, R0, t)
        vals = [f(P) - f(mn(R0)) for g in nonrh(rem) + nonrh(cre) for P in g[1]]
        lo = min(vals) if lo is None else min(lo, min(vals)); hi = max(vals) if hi is None else max(hi, max(vals))
amax = max(abs(f(a)) for a in ALPHA.values()); dmin = min(sum(p * d for p, d in zip(best, D[t])) for t in D)
N = (hi - lo + amax) // dmin + 1
print('window [%d,%d], |phi(alpha)|<=%d, step>=%d => tiles can coincide only if i-j <= %d' % (lo, hi, amax, dmin, N))
ok = True; cnt = 0
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL])
    for n in range(2, N + 2):
        for word in itertools.product(('a0', 'a1'), repeat=n):
            ws = walk(word, TE, R0)
            Rj, _, _ = ws[0]; Ri, _, Ci = ws[-1]
            Gj = set(nonrh(step(TE, tr(Rj, mal), word[0])[0]))
            cnt += 1
            if Gj & Ci: ok = False; print('HIT', TL, word)
print('(ii) words checked', cnt, 'no E^L_j tile in C(L_i), i>j:', ok)
