"""Merge' inductive step, absence form.  E^L sits at R - alpha (R = R_i, class 0), and E_i = L_i is the class-0 a-step t at R.
Known on E^L's board: nonrh C(L_i) (= nonrh C(E_i) by Lemma c0), created by J_L, final.
Claim: every type t' != t available to T_E at R - alpha has a removed tile overlapping (not equal to) a tile of nonrh C(E_i).
Then the board is a tiling => t' absent; E^L steps (hyp.) => by Lemma U it takes the present good type, which is t.
Also: G_t(R-alpha) compatible with C(E_i) (sanity), and which of its tiles ARE in C(E_i) (the 'CL+0' trail)."""
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/base.py').read().split("ALPHA = {}")[0])
ALPHA = {'B': (0, 0, 0, 1), 'C': (0, -1, 0, 0)}
ok_all = True
for TL in 'BC':
    TE = OTHER[TL]; al = ALPHA[TL]; mal = tuple(-a for a in al)
    R = canon(TE, 0); assert cls(TE, R) == 0 and canon(TL, 0) == R
    for t in ('a0', 'a1'):
        remT, RnT, creT = step(TE, R, t)
        remL, RnL, creL = step(TL, R, t); assert (set(remT), RnT, set(creT)) == (set(remL), RnL, set(creL)), 'c0'
        Ra = tr(R, mal); assert cls(TE, Ra) == 0
        K = set(nonrh(creT))
        cfg = configs_at(TE, Ra)
        Gt = set(nonrh(cfg[t][0]))
        line = []
        for tp, (rem, _, cre) in cfg.items():
            if tp == t: continue
            k = killers(set(nonrh(rem)), K)
            line.append((tp, k or 'NOT KILLED'))
            ok_all &= bool(k)
        print('T_L', TL, 'T_E', TE, 'step', t, '| wrong types:', line, '| G_t(R-a) compatible with C(E_i):', not killers(Gt, K),
              '| G_t(R-a) tiles in C(E_i):', sorted(g[0] for g in Gt & K), 'of', len(Gt),
              '| created rhombus of t at R-a == Rn - alpha:', cfg[t][1] == tr(RnT, mal))
print('ALL WRONG TYPES KILLED BY C(L_i):', ok_all)
