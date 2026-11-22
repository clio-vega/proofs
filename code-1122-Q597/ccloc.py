"""PROVE 11-22 Step 1b: (CC) as an exhaustive local check at the first contact.
Setting (all proved in the tex): at the divergence p, X_p = M0.prefix and Y_p = Mx.prefix differ exactly by W | Xadd
(Lemma Dinv); A takes the unique good T-configuration G_t present on X_p, the re-run takes the unique one G_t' present on
Y_p, t' != t (Lemma U); H_p meets W and G_t' meets Xadd (Theorem FC).
Unknowns enumerated (exhaustively, up to translation; R_p = canonical class-cA rhombus of T at base 0):
  cA, t, t'; the other journey's SEGMENT = its steps from the first to the last ANCHOR step, where an anchor step removes
  a tile of G_t (then that tile is in W) or creates a tile of G_t' (then it is in Xadd).  Every W tile of H_p and every
  Xadd tile of G_t' is anchored by EQUALITY with a removed/created tile of a good T'-configuration, so the anchor
  placements are finite; between first and last anchor the 4-offset om_T'(R_p, R^o_j) is monotone (Lemma window), so
  the segment length is bounded and the DFS below is exhaustive (it prunes only on the monotone box and on necessary
  conditions).  Optionally A's last h steps (backward words into R_p).
Necessary conditions imposed (each is a consequence of proved facts; tags used for the lesions):
  J1  the segment is a legal fragment: removed non-rh tiles were not created earlier in it (no revisit), created tiles
      were not removed earlier (Lemma nr), and the tiles known present on J_o's board at each time are pairwise
      equal-or-interior-disjoint;
  J2  (choice rule, J_o) at each segment step no OTHER good T'-configuration has all 3 tiles known present;
  XY  X_p contains R_p, G_t, all W tiles of the segment and the common part of G_t'; Y_p contains R_p, G_t' and all
      created tiles of the segment, and the common part of G_t: each family pairwise equal-or-interior-disjoint;
      S = G_t & W_seg != 0, V = G_t' & Xadd_seg != 0, G_t' & W_seg = 0, G_t & Xadd_seg = 0;
  UQ  (choice rule, A and re-run) no good T-configuration other than G_t (resp. G_t') has all tiles known on X_p (Y_p);
  HI  (A's history, h >= 1) the hexagons of A's last h steps lie in U_p, which is interior-disjoint from J_o's footprint
      (Lemma p2 proof): no segment tile (removed or created, rhombi included) meets them; and at each history step the
      present-set on X_k contains no extra certified configuration (UQ at k, W tiles are on X_k).
Output: every surviving (cA, t, t', segment); for each, c0 = first segment step removing a tile of G_t, and
  CC1: S <= rem(c0) | rem(c0+1);  CC2: (J class/type at c0, c0+1; A (cA,t); base(R_p) - base(R^o_c0)) in the 4-list.
usage: python3 ccloc.py H [LESION]   LESION in {none, noJ2, noUQ, noS2 (= noJ2+noUQ), noHI}"""
import sys, itertools
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
if os.environ.get('SAME', '0') == '1': OTHER = {'B': 'B', 'C': 'C'}   # control: both journeys to the same target (B/B)
H = int(sys.argv[1]); LES = sys.argv[2] if len(sys.argv) > 2 else 'none'
NOV = LES == 'noV'; USE_J2 = LES not in ('noJ2', 'noS2', 'noV'); USE_UQ = LES not in ('noUQ', 'noS2'); USE_HI = LES != 'noHI'
CCS = {('B', (0, 'a0'), (0, 'a1'), (1, 'a1'), (0, 0, 0, 0)), ('B', (1, 'b'), (0, 'a1'), (1, 'a1'), (0, 0, 1, -1)),
       ('C', (0, 'a1'), (0, 'a0'), (1, 'a0'), (0, 1, 0, 1)), ('C', (1, 'b'), (0, 'a0'), (1, 'a0'), (1, 0, 0, 0))}
st = Counter(); SURV = []
def shift_to(t_cfg, t_target):
    """translation v with tr(t_cfg, v) == t_target, or None (same kind and shape required)"""
    if t_cfg[0] != t_target[0]: return None
    v = tuple(a - b for a, b in zip(t_target[1][0], t_cfg[1][0]))
    return v if tr(t_cfg, v) == t_target else None
def placements(Tp, targets):
    """all (R^o) placements of a good T'-configuration whose removed non-rh tiles or created non-rh tiles contain a target"""
    out = set()
    for cJ in (0, 1):
        Rc = canon(Tp, cJ)
        for typ in TYPES:
            rem, Rn, cre = step(Tp, Rc, typ)
            for tc in nonrh(rem) + nonrh(cre):
                for tg in targets:
                    v = shift_to(tc, tg)
                    if v is not None: out.add(tr(Rc, v))
    return out
def preds(T, R):
    """(R_prev, removed, created) of every good T-step whose image is R"""
    out = []
    for c in (0, 1):
        Rc = canon(T, c)
        for typ in TYPES:
            rem, Rn, cre = step(T, Rc, typ)
            if classkey(Rn) != classkey(R): continue
            v = tuple(a - b for a, b in zip(base(R), base(Rn)))
            Rp = tr(Rc, v); rem2, Rn2, cre2 = step(T, Rp, typ)
            if Rn2 == R: out.append((Rp, typ, rem2, cre2))
    return out
def histories(T, R, h):
    if h == 0: yield []; return
    for Rp, typ, rem, cre in preds(T, R):
        for hist in histories(T, Rp, h - 1): yield [(Rp, typ, rem, cre)] + hist
for T in 'BC':
    Tp = OTHER[T]
    for cA in (0, 1):
        RA = canon(T, cA)
        G = {typ: nonrh(step(T, RA, typ)[0]) for typ in TYPES}
        allanch = placements(Tp, [x for typ in TYPES for x in G[typ]])
        Q = [offs(RA, R, Tp) for R in allanch]
        hi = [max(q[i] for q in Q) for i in (0, 1)]; lo = [min(q[i] for q in Q) for i in (2, 3)]
        D = max(sum(abs(a-b) for a, b in zip(q1, q2)) for q1 in Q for q2 in Q)
        st[(T, cA, 'anchor placements', len(allanch), 'l1-diam', D)] += 1
        HISTS = list(histories(T, RA, H)) if H > 0 else [[]]
        st[(T, cA, 'histories', len(HISTS))] += 1
        for t in TYPES:
            for t2 in TYPES:
                if t2 == t: continue
                Gt, Gt2 = G[t], G[t2]
                anch = placements(Tp, Gt + Gt2)
                maxdepth = [0]
                def inbox(R):
                    o = offs(RA, R, Tp)
                    return o[0] <= hi[0] and o[1] <= hi[1] and o[2] >= lo[0] and o[3] >= lo[1]
                def dfs(seg, R, known, created, removedW, removed_all):
                    """seg: list of (R_j, typ, rem, cre). known: tiles known on J_o's board now (W tiles later removed are
                    added when removed -- J1 is checked against everything known at the time of each removal)."""
                    maxdepth[0] = max(maxdepth[0], len(seg))
                    if seg:
                        yield seg, removedW, created
                    if not inbox(R): return
                    for typ in TYPES:
                        rem, Rn, cre = step(Tp, R, typ)
                        nr = nonrh(rem)
                        if any(x in created for x in nr): st['J1 prune: revisit'] += 1; continue
                        if any(x in removed_all for x in nonrh(cre)): st['J1 prune: re-creation'] += 1; continue
                        # W tiles removed now were on J_o's board throughout the segment: compatible with all known
                        newW = [x for x in nr]
                        if not all(compatible(x, list(known) + list(removedW)) for x in newW): st['J1 prune: overlap'] += 1; continue
                        o1, o2 = offs(RA, R, Tp), offs(RA, Rn, Tp)
                        mono = o2[0] >= o1[0] and o2[1] >= o1[1] and o2[2] <= o1[2] and o2[3] <= o1[3]
                        st[('ASSERT om monotone along J_o step (psi-part up, lambda-part down)', mono)] += 1
                        k2 = (set(known) - set(rem)) | set(nonrh(cre))
                        yield from dfs(seg + [(R, typ, rem, cre)], Rn, k2 | {Rn}, created | set(nonrh(cre)),
                                       removedW | set(nr), removed_all | set(rem))
                for R0 in anch:
                    # the segment's first step must be an anchor step: checked on output
                    for seg, Wseg, Xseg in dfs([], R0, {R0}, set(), set(), set()):
                        first, last = seg[0], seg[-1]
                        isanch = lambda s: bool(set(nonrh(s[2])) & set(Gt)) or bool(set(nonrh(s[3])) & set(Gt2))
                        if not (isanch(first) and isanch(last)): continue
                        if USE_J2:
                            # choice rule for J_o, post hoc with full segment knowledge: at step j its board contains R_j,
                            # the tiles created at steps < j (never removed later: no revisit), and every tile removed at a
                            # step >= j (an M0 tile not yet removed).  No other good T'-configuration may be all there.
                            badj = False
                            for j, (Rj, typj, remj, crej) in enumerate(seg):
                                known_j = {Rj} | {x for s_ in seg[:j] for x in nonrh(s_[3])} | {x for s_ in seg[j:] for x in nonrh(s_[2])}
                                for o in TYPES:
                                    if o == typj: continue
                                    if all(x in known_j for x in nonrh(step(Tp, Rj, o)[0])): badj = True
                            if badj: st['J2 prune'] += 1; continue
                        S = set(Gt) & Wseg; V = set(Gt2) & Xseg
                        if NOV: V = V or set()
                        if not S or (not V and not NOV): st['XY: S or V empty'] += 1; continue
                        if (set(Gt2) & Wseg and not NOV) or set(Gt) & Xseg: st['XY: Gt2 meets W or Gt meets Xadd'] += 1; continue
                        comG = [x for x in Gt2 if x not in V]; comH = [x for x in Gt if x not in S]
                        if NOV: comG = []
                        Xk = [RA] + Gt + sorted(Wseg) + comG; Yk = [RA] + (Gt2 if not NOV else []) + sorted(Xseg) + comH
                        okXY = all(a == b or not overlap(a, b) for a, b in itertools.combinations(Xk, 2)) and \
                               all(a == b or not overlap(a, b) for a, b in itertools.combinations(Yk, 2))
                        if not okXY: st['XY prune'] += 1; continue
                        if USE_UQ:
                            if any(all(x in Xk for x in G[o]) for o in TYPES if o != t): st['UQ prune X'] += 1; continue
                            if not NOV and any(all(x in Yk for x in G[o]) for o in TYPES if o != t2): st['UQ prune Y'] += 1; continue
                        foot = [x for s in seg for x in s[2]] + [x for s in seg for x in s[3]]
                        survivors_h = 0
                        for hist in HISTS:
                            if USE_HI and hist:
                                hexes = [x for (_, _, rem, cre) in hist for x in rem]
                                if any(meets(x, hexes) for x in foot): continue
                                # UQ at history steps on X_k: W tiles and the history step's own hexagon are on X_k
                                okh = True
                                for (Rk, typk, remk, crek) in hist:
                                    for o in TYPES:
                                        if o == typk: continue
                                        Go = nonrh(step(T, Rk, o)[0])
                                        if all(x in Wseg or x in remk for x in Go): okh = False
                                if not okh: continue
                            survivors_h += 1
                        if survivors_h == 0: st['HI prune'] += 1; continue
                        # outcome
                        c0 = next(j for j, s in enumerate(seg) if set(nonrh(s[2])) & set(Gt))
                        rd = lambda j: set(nonrh(seg[j][2])) if j < len(seg) else set()
                        cc1 = S <= rd(c0) | rd(c0 + 1)
                        tj = lambda j: (cls(Tp, seg[j][0]), seg[j][1]) if j < len(seg) else None
                        off = tuple(a - b for a, b in zip(base(RA), base(seg[c0][0])))
                        cc2 = (T, tj(c0), tj(c0 + 1), (cA, t), off) in CCS
                        sheet_p = not any(set(nonrh(s_[2])) & S for s_ in seg if cls(Tp, s_[0]) == 1 and s_[1] in ('a0', 'a1'))
                        st[('EQUIV', T, 'Sheet_p', sheet_p, 'CC', cc1 and cc2)] += 1
                        key = (T, 'A', (cA, t), "re-run", t2 if not NOV else '-', 'Sheet_p', sheet_p, 'J c0,c0+1', tj(c0), tj(c0 + 1), 'off', off, 'CC1', cc1, 'CC2', cc2,
                               '|S|', len(S), '|V|', len(V))
                        st[('SURV',) + key] += 1
                        if st[('SURV',) + key] == 1: SURV.append((key, [(cls(Tp, s[0]), s[1]) for s in seg]))
                st[(T, cA, t, t2, 'max segment depth')] = max(st[(T, cA, t, t2, 'max segment depth')], maxdepth[0])
print('H', H, 'LESION', LES)
for k in sorted(st, key=str): print(st[k], k)
print('examples:')
for key, w in SURV: print(key, 'segment', w)
