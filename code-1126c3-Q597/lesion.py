"""Lesion of stepchk: replace alpha by each other unit vector v of Z^4.  Does 'every wrong type killed by C(L_i) AND
G_t(R-v) compatible with C(L_i)' still hold?  If the step lemma passes for many v it is not using alpha."""
import itertools
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/stepchk.py').read().split("ok_all = True")[0])
UNITS = [tuple(int(i == j) * s for j in range(4)) for i in range(4) for s in (1, -1)]
for TL in 'BC':
    TE = OTHER[TL]
    for v in UNITS:
        mal = tuple(-a for a in v); R = canon(TE, 0); good = True
        for t in ('a0', 'a1'):
            remT, RnT, creT = step(TE, R, t); K = set(nonrh(creT)); Ra = tr(R, mal)
            if cls(TE, Ra) != 0: good = False; break
            cfg = configs_at(TE, Ra)
            if killers(set(nonrh(cfg[t][0])), K): good = False
            for tp, (rem, _, _) in cfg.items():
                if tp != t and not killers(set(nonrh(rem)), K): good = False
        print('T_L', TL, 'v', v, 'is alpha' if v == ALPHA[TL] else '', 'step lemma holds:', good)
