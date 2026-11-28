"""Frozen interpolation.  E_T(j): pop x, take j steps of x toward T, then run y's whole journey toward the OTHER nest,
then resume x toward T.  E_B(len K) = L, E_C(len H) = R.  Is E_T(j) independent of j (where defined)?  Is E_B(j)=E_C(j')?
usage: python3 frozen.py NSEED"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
def packed(b, R):
    vs = b.tiles[R][1]
    for X, info in enumerate(b.nest):
        if X == 0: continue
        ij = b.packed_at(info, vs)
        if ij is not None and ij not in info['pos'] and b.is_ideal(set(info['pos']) | {ij}): return True
    return False
def E(b, x, y, dx, dy, j):
    b = copy.deepcopy(b); R = b.nest[0]['pos'].pop(x); d = unit_dir(dx)
    for _ in range(j):
        if packed(b, R): return 'x landed before j'
        R = b.step(R, d)
        if R is None: return 'x stuck'
    try:
        t = b.nest[0]['pos'].pop(y); b.journey(t, 0, unit_dir(dy))
    except RuntimeError as e: return 'y ' + str(e)[:20]
    try: b.journey(R, 0, d)
    except RuntimeError as e: return 'x resume ' + str(e)[:20]
    return tuple(tiles_of(b))
st = Counter(); NSEED = int(sys.argv[1])
for k, (key, b0) in enumerate(census()):
    if k >= NSEED: break
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S)
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            st['SW n'] += 1
            res = {}
            for T, O in (('B', 'C'), ('C', 'B')):
                res[T] = [E(b, x, y, dirs[T], dirs[O], j) for j in range(0, 12)]
            Lt = res['B'][-1] if False else None
            # full values: find last j where tuple result
            for T in 'BC':
                vals = [v for v in res[T] if not isinstance(v, str)]
                errs = Counter(v for v in res[T] if isinstance(v, str))
                for e in errs: st[(T, 'err', e)] += 1
                st[(T, 'distinct results over j')] += 0
                st[(T, 'n distinct', len(set(vals)))] += 1
                # which j agree with the final (landed) value
                fin = vals[-1]
                pat = ''.join('.' if isinstance(v, str) else ('=' if v == fin else 'x') for v in res[T])
                st[(T, 'pattern', pat)] += 1
            st['E_B final == E_C final (SW)'] += [v for v in res['B'] if not isinstance(v, str)][-1] == [v for v in res['C'] if not isinstance(v, str)][-1]
for k in sorted(st, key=str): print(st[k], k)
