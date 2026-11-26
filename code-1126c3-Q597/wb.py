"""(W_B) first step: E^L at R_{k+2}-alpha (class 0, no stop).  Known present: nonrh C(L_{k+2}), L_{k+2} = T_L's (0,b) at R_{k+2}.
Which T_E types at R-alpha survive (no removed tile overlapping a known tile)?  Also with K += nonrh C(E^L_{k+1}) for each a-type
of the last corridor step (just created by E^L, certainly present)."""
exec(open('/home/clio/projects/proofs/code-1126c3-Q597/stepchk.py').read().split("ok_all = True")[0])
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remLb, RnLb, creLb = step(TL, R, 'b'); K = set(nonrh(creLb))
    cfg = configs_at(TE, Ra)
    surv = [tp for tp, (rem, _, _) in cfg.items() if not killers(set(nonrh(rem)), K)]
    print('T_L', TL, 'survivors with C(L_{k+2}) only:', surv)
    # previous corridor step t0 creating R: find R0 with step(TE,R0,t0) creating R
    for t0 in ('a0', 'a1'):
        # invert: candidate predecessor = translate so that created rhombus equals R
        rem0, Rn0, cre0 = step(TE, R, t0); d = tuple(a - b for a, b in zip(min(R[1]), min(Rn0[1])))
        Rp = tr(R, d); rem0, Rn0, cre0 = step(TE, Rp, t0); assert Rn0 == R
        K2 = K | {tr(g, mal) for g in nonrh(cre0)}
        surv2 = [tp for tp, (rem, _, _) in cfg.items() if not killers(set(nonrh(rem)), K2)]
        print('    prev corridor step', t0, 'survivors:', surv2)
    # what do the survivors look like vs E_{k+2} = T_E (0,b) at R translated
    remEb, RnEb, creEb = step(TE, R, 'b')
    for tp in surv:
        print('    survivor', tp, 'is E_{k+2}-alpha:', step(TE, Ra, tp)[1] == tr(RnEb, mal) and tp == 'b',
              '| created rhombus class', cls(TE, step(TE, Ra, tp)[1]))
print('---- E^L_{k+3}')
aT = {'B': 'a1', 'C': 'a0'}
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remLb, RnLb, creLb = step(TL, R, 'b'); K = set(nonrh(creLb))
    rem2, R3, cre2 = step(TE, Ra, aT[TL])
    K3 = (K - set(nonrh(rem2))) | set(nonrh(cre2))
    cfg = configs_at(TE, R3)
    surv = [tp for tp, (rem, _, _) in cfg.items() if not killers(set(nonrh(rem)), K3)]
    print('T_L', TL, 'R3 class', cls(TE, R3), 'survivors:', surv,
          {tp: sorted(('C(Lk+2)' if g in K else 'C(EL)' if g in set(cre2) else 'other') for g in nonrh(cfg[tp][0])) for tp in surv})
print('---- coupling E_{k+3} (J_E, on M0) vs E^L_{k+3} candidates')
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remEb, RE3, creEb = step(TE, R, 'b')
    rem2, R3, cre2 = step(TE, Ra, aT[TL])
    cE = configs_at(TE, RE3); cEL = configs_at(TE, R3)
    for x, (rx, _, _) in cE.items():
        Gx = set(nonrh(rx))
        for y in (aT[TE], 'b'):
            Gy = set(nonrh(cEL[y][0]))
            print('T_L', TL, 'E_{k+3}=(%d,%s)' % (cls(TE, RE3), x), 'E^L_{k+3}=(0,%s)' % y, 'equal tiles', sorted(g[0] for g in Gx & Gy),
                  'overlapping pairs', killers(Gx, Gy), '| E_{k+3} tiles vs C(E_{k+2}):', len(Gx & set(creEb)))
print('---- bad pairs killed by an overlap avoiding the C(L_{k+2}) tile of E^L_{k+3}?')
allok = True
for TL in 'BC':
    TE = OTHER[TL]; mal = tuple(-a for a in ALPHA[TL]); R = canon(TE, 0); Ra = tr(R, mal)
    remLb, RnLb, creLb = step(TL, R, 'b'); CL = set(creLb)
    remEb, RE3, creEb = step(TE, R, 'b'); rem2, R3, cre2 = step(TE, Ra, aT[TL])
    cE = configs_at(TE, RE3); cEL = configs_at(TE, R3)
    for x, (rx, _, _) in cE.items():
        Gx = set(nonrh(rx))
        for y in (aT[TE], 'b'):
            Gy = set(nonrh(cEL[y][0])) - CL
            good = (x, y) in {('b', aT[TE]), (aT[TL], 'b')}
            k = killers(Gx, Gy)
            print('T_L', TL, x, y, 'census pair' if good else 'bad pair', 'M0-only overlaps:', k or 'none')
            allok &= (not k) if good else bool(k)
print('PAIR LEMMA (M0 tiles only) OK:', allok)
