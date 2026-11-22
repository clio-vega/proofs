"""PROVE 11-22 Step 1c: the local contact check with PADDING.
Core: ccloc.py's anchor segments (exhaustive).  Each core survivor (T, cA, t, t', seg) is extended by
  * KB steps of J_o BEFORE the segment (all legal predecessor words), or fewer if J_o STARTS there (start rhombus class 0;
    it is then an M0 tile in W -- the nest-A rhombus of J_o's corner);
  * KA steps AFTER the segment, or fewer if J_o ENDS there (final rhombus class 1, Lemma fin; it is then in Xadd);
  * A's last H steps (all backward words into R_p; H = 0 allowed).
Padding steps are NON-anchor (they neither remove a tile of G_t nor create a tile of G_t'): the segment runs from the
first to the last anchor by definition.  Every real contact therefore appears as (core, some padding, some history),
truncated, so a configuration killed for EVERY padding and history is not realisable.
Necessary conditions checked on the extended picture (all from proved facts):
  J1  legality of J_o's word: no revisit, no re-creation; at each time the tiles known on J_o's board (R_j, earlier
      created tiles, tiles removed at >= j) are pairwise equal-or-interior-disjoint.
  J2  choice rule for J_o at every step (no other good T'-configuration fully known present).
  XY  X_p contains R_p, G_t, every W tile of the extended word (incl. a start rhombus), the common part of G_t', and A's
      created tiles still present (history); Y_p contains R_p, G_t', every created tile of the extended word (incl. a final
      rhombus), the common part of G_t, and A's present created tiles: each family pairwise equal-or-interior-disjoint.
  UQ  choice rule for A on X_p and for the re-run on Y_p.
  HI  A's history hexagons are interior-disjoint from J_o's whole extended footprint (U_p vs F_o, Lemma p2 proof), and
      the choice rule for A at each history step on what is known of X_k (W tiles, A's own tiles present then).
usage: python3 ccpad.py KB KA H [LESION]  (LESION: none | noS2 = drop J2 and UQ everywhere | noPAD = KB=KA=0)"""
import sys, itertools
sys.argv = [sys.argv[0], '0', 'none'] + sys.argv[1:]
KB, KA, H = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]); LES2 = sys.argv[6] if len(sys.argv) > 6 else 'none'
_src = open('/home/clio/projects/proofs/code-1122-Q597/ccloc.py').read()
_src = _src.replace("print('H', H, 'LESION', LES)", "CORE = SURV_ALL\nif False:\n    pass").replace("for k in sorted(st, key=str): print(st[k], k)\nprint('examples:')\nfor key, w in SURV: print(key, 'segment', w)", "")
_src = _src.replace("                        if st[('SURV',) + key] == 1: SURV.append", "                        SURV_ALL.append((T, cA, t, t2, list(seg), key))\n                        if st[('SURV',) + key] == 1: SURV.append")
_src = _src.replace("st = Counter(); SURV = []", "st = Counter(); SURV = []; SURV_ALL = []")
if LES2 == 'noS2': _src = _src.replace("LES = sys.argv[2] if len(sys.argv) > 2 else 'none'", "LES = 'noS2'")
exec(_src)
KB, KA, H = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
USE_S2 = LES2 != 'noS2'
NOSE = os.environ.get('NOSE', '0') == '1'   # diagnostic: forbid J_o start/end inside the padding
st2 = Counter()
def pw(a, b): return all(x == y or not overlap(x, y) for x, y in itertools.combinations(list(a) + list(b), 2))
def jpreds(Tp, R):
    out = []
    for c in (0, 1):
        Rc = canon(Tp, c)
        for typ in TYPES:
            rem, Rn, cre = step(Tp, Rc, typ)
            if classkey(Rn) != classkey(R): continue
            v = tuple(a - b for a, b in zip(base(R), base(Rn)))
            Rp_ = tr(Rc, v); rem2, Rn2, cre2 = step(Tp, Rp_, typ)
            if Rn2 == R: out.append((Rp_, typ, rem2, cre2))
    return out
def backs(Tp, R, k):
    """backward extensions of length k (list of steps, earliest first) or shorter ones ending at a START (class 0)"""
    if cls(Tp, R) == 0 and not NOSE: yield ([], True)
    if k == 0: yield ([], False); return
    for s in jpreds(Tp, R):
        for ext, startf in backs(Tp, s[0], k - 1): yield (ext + [s], startf)
def fwds(Tp, R, k):
    if cls(Tp, R) == 1 and not NOSE: yield ([], True)
    if k == 0: yield ([], False); return
    for typ in TYPES:
        rem, Rn, cre = step(Tp, R, typ); s = (R, typ, rem, cre)
        for ext, endf in fwds(Tp, Rn, k - 1): yield ([s] + ext, endf)
def full_check(T, cA, t, t2, word, startf, endf, hist):
    Tp = OTHER[T]; RA = canon(T, cA)
    G = {typ: nonrh(step(T, RA, typ)[0]) for typ in TYPES}; Gt, Gt2 = G[t], G[t2]
    # J1 on the extended word
    created, removed_all = set(), set()
    for j, (Rj, typj, remj, crej) in enumerate(word):
        if any(x in created for x in nonrh(remj)): return 'J1 revisit'
        if any(x in removed_all for x in nonrh(crej)): return 'J1 re-creation'
        created |= set(nonrh(crej)); removed_all |= set(remj)
    Wext = {x for s in word for x in nonrh(s[2])} | ({word[0][0]} if startf else set())
    Xext = {x for s in word for x in nonrh(s[3])} | ({word[-1][3][-1]} if endf else set())
    for j, (Rj, typj, remj, crej) in enumerate(word):
        known_j = {Rj} | {x for s in word[:j] for x in nonrh(s[3])} | {x for s in word[j:] for x in nonrh(s[2])}
        if startf: known_j |= set()
        if not pw(known_j, []): return 'J1 overlap'
        if USE_S2:
            for o in TYPES:
                if o != typj and all(x in known_j for x in nonrh(step(Tp, Rj, o)[0])): return 'J2'
    S = set(Gt) & Wext; V = set(Gt2) & Xext
    if not S or not V: return 'XY S/V empty'
    if set(Gt2) & Wext or set(Gt) & Xext: return 'XY cross'
    # A's history: present created tiles at time p, and hexagons
    Apres = set(); hexes = []
    for k, (Rk, typk, remk, crek) in enumerate(hist):       # hist[0] is step p-1, hist[1] step p-2, ...
        hexes += remk
    # tiles present at time p from A's last H steps: created at p-1; created at p-2 unless removed at p-1; ...
    for k in range(len(hist)):
        later_removed = {x for kk in range(k) for x in hist[kk][2]}
        Apres |= {x for x in nonrh(hist[k][3]) if x not in later_removed}
    comG = [x for x in Gt2 if x not in V]; comH = [x for x in Gt if x not in S]
    Xk = {RA} | set(Gt) | Wext | set(comG) | Apres; Yk = {RA} | set(Gt2) | Xext | set(comH) | Apres
    if not pw(Xk, []): return 'XY X overlap'
    if not pw(Yk, []): return 'XY Y overlap'
    if USE_S2:
        if any(all(x in Xk for x in G[o]) for o in TYPES if o != t): return 'UQ X'
        if any(all(x in Yk for x in G[o]) for o in TYPES if o != t2): return 'UQ Y'
    if hist:
        foot = [x for s in word for x in s[2]] + [x for s in word for x in s[3]]
        if any(meets(x, hexes) for x in foot): return 'HI footprint meets U_p'
        if USE_S2:
            # choice rule for A at history step p-1-k on X_{p-1-k}: W tiles are there, and that step's own hexagon
            for k, (Rk, typk, remk, crek) in enumerate(hist):
                for o in TYPES:
                    if o != typk and all(x in Wext or x in remk for x in nonrh(step(T, Rk, o)[0])): return 'HI UQ'
    return 'OK'
CCSET = CCS
res = Counter(); EX = {}
for (T, cA, t, t2, seg, key) in CORE:
    Tp = OTHER[T]; RA = canon(T, cA)
    G = {typ: nonrh(step(T, RA, typ)[0]) for typ in TYPES}
    isanch = lambda s: bool(set(nonrh(s[2])) & set(G[t])) or bool(set(nonrh(s[3])) & set(G[t2]))
    BK = [(e, f) for e, f in backs(Tp, seg[0][0], KB) if not any(isanch(s) for s in e)]
    FW = [(e, f) for e, f in fwds(Tp, seg[-1][2] and step(Tp, seg[-1][0], seg[-1][1])[1], KA) if not any(isanch(s) for s in e)]
    HS = list(histories(T, RA, H)) if H > 0 else [[]]
    why = Counter(); ok = 0
    for (be, sf) in BK:
        for (fe, ef) in FW:
            word = be + seg + fe
            for hist in HS:
                r = full_check(T, cA, t, t2, word, sf, ef, hist)
                why[r] += 1
                if r == 'OK':
                    ok += 1
                    if (key, 'ex') not in EX: EX[(key, 'ex')] = ([(cls(Tp, s[0]), s[1]) for s in word], sf, ef, [(cls(T, h[0]), h[1]) for h in hist])
    res[(key, 'realisable (some padding+history passes)', ok > 0)] += 1
    for r, n in why.items(): st2[(key[:10], r)] += n
print('KB', KB, 'KA', KA, 'H', H, 'LESION', LES2)
for k in sorted(res, key=str): print(res[k], k)
print('--- survivors with an example extended word (J_o word, startflag, endflag, A history newest-first)')
for k, v in EX.items(): print(k[0][:12], v)
print('--- reasons (per core key, summed over paddings x histories)')
for k in sorted(st2, key=str): print(st2[k], k)
