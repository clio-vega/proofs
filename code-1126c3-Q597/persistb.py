"""Persistence for the exit: tiles created by L_{k+2} = T_L's (0,b) at R_{k+2} are not removed by E^L_j = E_j - alpha
(2<=j<=k+1, class-0 a-steps), nor by E^L_0, E^L_1 (those remove C(L_0),C(L_1),N_E tiles: created-once / M0 argument).
Window as in persist.py, with the created tiles of the b-step included in the window."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/persist.py').read().split("ok = True; cnt = 0")[0])
lo2, hi2 = lo, hi
for T in 'BC':
    rem, Rn, cre = step(T, R0, 'b')
    vals = [f(P) - f(mn(R0)) for g in nonrh(cre) for P in g[1]]
    lo2, hi2 = min(lo2, min(vals)), max(hi2, max(vals))
N2 = (hi2 - lo2 + amax) // dmin + 1
print('window with b-created tiles [%d,%d] => coincidence needs i-j <= %d' % (lo2, hi2, N2))
ok = True; cnt = 0
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL])
    for n in range(1, N2 + 1):
        for word in itertools.product(('a0', 'a1'), repeat=n):
            ws = walk(word, TE, R0)
            Rj = ws[0][0]
            R = R0
            for t in word: R = step(TE, R, t)[1]
            Cb = set(nonrh(step(TL, R, 'b')[2]))
            Gj = set(nonrh(step(TE, tr(Rj, mal), word[0])[0])); cnt += 1
            if Gj & Cb: ok = False; print('HIT', TL, word)
print('words', cnt, 'no E^L_j tile in C(L_{k+2}):', ok)
