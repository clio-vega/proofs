"""Lemma sq without the segment enumeration?  At the contact rhombus R (class cA of T), A takes G_t, the re-run G_t' (t != t').
Unknowns: S = G_t cap W (nonempty), V = G_t' cap Xadd (nonempty) -- Theorem fc.  Proved constraints (no J_o segment):
  (i)   S cap G_t' = {} and V cap G_t = {}   (W tiles are not on Y_p, Xadd tiles are not on X_p; Lemma dinv);
  (ii)  X_p contains R, G_t and G_t' minus V; Y_p contains R, G_t' and G_t minus S: each family pairwise equal-or-disjoint;
  (iii) Lemma U on X_p and on Y_p: no third type has all three tiles in the known part.
Question: is there a solution with the square of G_t NOT in S?  Also prints which tiles G_t and G_t' share."""
import itertools
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
st = Counter()
def okfam(F): return all(a == b or not overlap(a, b) for a, b in itertools.combinations(F, 2))
def subsets(L): return [set(c) for r in range(1, len(L) + 1) for c in itertools.combinations(L, r)]
for T in 'BC':
    for cA in (0, 1):
        R = canon(T, cA); G = {t: nonrh(step(T, R, t)[0]) for t in TYPES}
        for t in TYPES:
            for t2 in TYPES:
                if t2 == t: continue
                Gt, Gt2 = G[t], G[t2]
                shared = sorted(u[0] for u in set(Gt) & set(Gt2))
                sq = next(u for u in Gt if u[0] == 'sq')
                sols = []
                for S in subsets(Gt):
                    if S & set(Gt2): continue
                    for V in subsets(Gt2):
                        if V & set(Gt): continue
                        X = [R] + Gt + [u for u in Gt2 if u not in V]; Y = [R] + Gt2 + [u for u in Gt if u not in S]
                        if not (okfam(X) and okfam(Y)): continue
                        if any(all(g in X for g in G[o]) for o in TYPES if o != t): continue
                        if any(all(g in Y for g in G[o]) for o in TYPES if o != t2): continue
                        sols.append((tuple(sorted(u[0] for u in S)), tuple(sorted(u[0] for u in V)), sq in S))
                nosq = [s for s in sols if not s[2]]
                st[(T, cA, t, t2, 'shared', tuple(shared), 'solutions', len(sols), 'with sq notin S', len(nosq), tuple(sorted(set(nosq))))] += 1
for k in sorted(st, key=str): print(st[k], k)
# The by-hand route: for t != t', the square of G_t is not a tile of G_t' and overlaps (interior) some tile of G_t'.
# Control: the triangles of G_t do NOT all have this property (else S could be anything).
print('--- overlap fact')
for T in 'BC':
    for cA in (0, 1):
        R = canon(T, cA); G = {t: nonrh(step(T, R, t)[0]) for t in TYPES}
        for t in TYPES:
            for t2 in TYPES:
                if t2 == t: continue
                sq = next(u for u in G[t] if u[0] == 'sq')
                hit = [g[0] for g in G[t2] if g != sq and overlap(sq, g)]
                tri = [(any(g != u and overlap(u, g) for g in G[t2]), u in G[t2]) for u in G[t] if u[0] == 'tri']
                print('OVL', T, cA, t, t2, 'sq in G_t2', sq in G[t2], 'sq overlaps', hit, 'tri (overlaps, shared)', tri)
