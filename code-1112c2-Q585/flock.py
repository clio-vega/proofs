"""Flocks, LR tableaux in a frame, and Purbhoo's Cor 3.4 swap."""
from migrate import *
import math, cmath

# ---------------- LR tableaux on a box set (French: (col,row), row up)
def lr_tableaux(boxes, content):
    boxes = sorted(boxes, key=lambda b: (b[1], -b[0]))     # rows south->north, each east->west = reversed reading word
    out = []; T = {}; cnt = [0]*(len(content)+2); need = list(content)
    def rec(k):
        if k == len(boxes):
            if cnt[1:len(content)+1] == need: out.append(dict(T))
            return
        c, r = boxes[k]
        lo = 1; hi = len(content)
        if (c, r-1) in T: lo = max(lo, T[(c, r-1)]+1)
        if (c+1, r) in T: hi = min(hi, T[(c+1, r)])
        for v in range(lo, hi+1):
            if cnt[v] >= need[v-1]: continue
            if v > 1 and cnt[v-1] <= cnt[v]: continue
            T[(c, r)] = v; cnt[v] += 1; rec(k+1); cnt[v] -= 1; del T[(c, r)]
    rec(0); return out
def is_lr(T):
    for (c, r), v in T.items():
        if (c+1, r) in T and T[(c+1, r)] < v: return False
        if (c, r+1) in T and T[(c, r+1)] <= v: return False
    rows = sorted(set(r for _, r in T), reverse=True)
    word = [T[(c, r)] for r in rows for c in sorted(c for c, rr in T if rr == r)]
    cnt = Counter()
    for v in reversed(word):
        cnt[v] += 1
        if v > 1 and cnt[v] > cnt[v-1]: return False
    return True
def content(T):
    c = Counter(T.values()); return tuple(c[k] for k in range(1, max(c)+1)) if c else ()
def skew_of(boxes):
    """normalise translation; return (outer, inner) partitions (row 0 = south) or None if not a skew shape"""
    if not boxes: return ((), ())
    c0 = min(c for c, r in boxes); r0 = min(r for c, r in boxes)
    B = {(c-c0, r-r0) for c, r in boxes}; R = max(r for c, r in B)+1
    outer = []; inner = []
    for r in range(R):
        cs = sorted(c for c, rr in B if rr == r)
        if not cs: return None
        if cs != list(range(cs[0], cs[-1]+1)): return None
        outer.append(cs[-1]+1); inner.append(cs[0])
    if any(outer[i] < outer[i+1] or inner[i] < inner[i+1] for i in range(R-1)): return None
    return tuple(outer), tuple(inner)

# ---------------- frames
def unit_dir(v): z = G.phi(v); return z/abs(z)
def box_of(board, X, vs, E, N):
    info = board.nest[X]
    c = tuple(sum(Fr(P[k]) for P in vs)/4 for k in range(4)); D = psub(c, info['X'])
    DE = tuple(D[k] if E[k] != 0 else 0 for k in range(4)); DN = psub(D, DE)
    x = Board.coef(DE, E) if any(DE) else Fr(0); y = Board.coef(DN, N) if any(DN) else Fr(0)
    assert x is not None and y is not None and x.denominator == 2 and y.denominator == 2, (x, y)
    return (math.floor(x), math.floor(y))
def induced(board, X, E, N):
    """orientation of nest X obtained by rotating (E,N) by <= 60 degrees"""
    u, w = board.nest[X]['u'], board.nest[X]['w']; neg = lambda v: tuple(-a for a in v)
    res = []
    for E2, N2 in ((u, w), (neg(u), neg(w)), (w, u), (neg(w), neg(u))):
        a1 = cmath.phase(unit_dir(E2)/unit_dir(E)); a2 = cmath.phase(unit_dir(N2)/unit_dir(N))
        if abs(a1-a2) < 1e-9 and abs(a1) <= math.pi/3+1e-9: res.append((E2, N2, a1))
    assert len(res) == 1, res
    return res[0]

def migrate_flock(board, src, flock, E, N, direction):
    """flock: dict tid->value (tiles in nest src); (E,N) its orientation; direction one of 'N','E','S','W'.
    Returns (dst, new flock dict tid->value, (E',N'), boxes dict box->value in induced frame, step counter)."""
    dvec = {'E': E, 'N': N, 'W': tuple(-a for a in E), 'S': tuple(-a for a in N)}[direction]
    d = unit_dir(dvec)
    key = {t: (v, box_of(board, src, board.tiles[t][1], E, N)[0]) for t, v in flock.items()}
    order = sorted(flock, key=lambda t: key[t], reverse=(direction in 'NE'))
    inv = {t: ij for ij, t in board.nest[src]['pos'].items()}
    out = {}; dst = None; tot = Counter()
    for t in order:
        del board.nest[src]['pos'][inv[t]]
        assert Board.is_ideal(set(board.nest[src]['pos'])) or True
        R, X, ij, st = board.journey(t, src, d); tot += st
        assert dst is None or dst == X, 'Prop 3.1 violated'
        dst = X; out[R] = flock[t]
    E2, N2, ang = induced(board, dst, E, N)
    boxes = {box_of(board, dst, board.tiles[R][1], E2, N2): v for R, v in out.items()}
    return dst, out, (E2, N2), boxes, tot

# ---------------- canonical mosaic (∅, γ∨, γ) with nest B = lam in frame (E_B,N_B)
def all_strings(n, k1):
    for S in itertools.combinations(range(n), k1):
        yield tuple(1 if i in S else 0 for i in range(n))
def nest_shape_of_string(n, side, string):
    """build only the nest rhombi of one side from a string (via a fake lift along that side)"""
    left, right, bottom = side_edges(n); edges = (left, right, bottom)[side]
    vec = (edges[0][1][0]-edges[0][0][0], edges[0][1][1]-edges[0][0][1])
    step = lambda l: (vec[0], vec[1], 0, 0) if l == 0 else (0, 0, vec[0], vec[1])
    start = (0, 0, 0, 0); target = sorted(string); cur = list(string); rh = []
    changed = True
    while changed:
        changed = False
        for k in range(n-1):
            if cur[k] != target[0] and cur[k+1] == target[0]:
                p = start
                for l in cur[:k]: p = padd(p, step(l))
                x, y = cur[k], cur[k+1]
                rh.append((p, padd(p, step(x)), padd(padd(p, step(x)), step(y)), padd(p, step(y))))
                cur[k], cur[k+1] = y, x; changed = True
    cnt = sum(1 for l in target if l == target[0]); X = start
    for l in target[:cnt]: X = padd(X, step(l))
    info = dict(X=X, u=tuple(-a for a in step(target[0])), w=step(target[-1]))
    fake = Board.__new__(Board)
    return {fake.packed_at(info, vs) for vs in rh}
def canonical(lam, k0, k1):
    n = k0+k1; left, right, bottom = side_edges(n)
    want = {(i, j) for j, l in enumerate(lam) for i in range(l)}
    strs = [s for s in all_strings(n, k1) if nest_shape_of_string(n, 1, s) == want]
    assert len(strs) == 1, len(strs)
    fixed = {}
    for e, l in zip(left, sorted([0]*k0+[1]*k1)): fixed[E(*e)] = l
    for e, l in zip(right, strs[0]): fixed[E(*e)] = l
    Ps = puzzles(n, fixed, limit=2)
    assert len(Ps) == 1, ('Lemma 2.6 uniqueness', len(Ps))
    M = build_mosaic(n, Ps[0], SIGN, RULE)
    b = Board(M)
    assert not b.nest[0]['pos'] and set(b.nest[1]['pos']) == want
    return b

def cor34(f, lam, k0, k1, verbose=False):
    """f: dict (col,row)->value, LR tableau on lam/mu (French, row 0 = longest).  Returns G' tableau (boxes in induced frame)."""
    b = canonical(lam, k0, k1); B = b.nest[1]
    EB, NB = B['u'], B['w']; neg = lambda v: tuple(-a for a in v)
    F = {B['pos'][(c, r)]: v for (c, r), v in f.items()}
    mu_pos = [ij for ij in B['pos'] if ij not in f]
    # G: unique flock on the remainder, orientation (-E_B,-N_B)
    boxesG = {box_of(b, 1, b.tiles[B['pos'][ij]][1], neg(EB), neg(NB)): ij for ij in mu_pos}
    mu_shape = tuple(sorted(Counter(r for _, r in mu_pos).values(), reverse=True))
    Ls = lr_tableaux(set(boxesG), mu_shape) if mu_pos else [{}]
    assert len(Ls) == 1, ('Lemma 2.5 uniqueness', len(Ls))
    Gfl = {B['pos'][boxesG[bx]]: v for bx, v in Ls[0].items()}
    d1, F2, fr1, boxF, st1 = migrate_flock(b, 1, F, EB, NB, 'E')
    assert is_lr(boxF), 'Thm 1 fails for F'
    if not Gfl: return {}, (d1, None), b
    d2, G2, fr2, boxG, st2 = migrate_flock(b, 1, Gfl, neg(EB), neg(NB), 'W')
    assert is_lr(boxG), 'Thm 1 fails for G'
    if verbose: print('F ->', d1, dict(st1), ' G ->', d2, dict(st2))
    return boxG, (d1, d2, boxF), b
