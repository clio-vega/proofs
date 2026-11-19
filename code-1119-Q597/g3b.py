"""PROVE 11-19 retreat G3 test (b), <=5. For the OTHER (non-clean) journey of each nontrivial B/C record:
is the pair of middles a two-block swap at some cut s: mB = (mA[s:] + v1) . (mA[:s] + v2), each block a rigid translate?
s=1 is forward BS (known 0/521); s=len-1 is BS of the inverted (reversed) journey.  Control: the clean journey (expect s=1)."""
import sys
exec(open('/home/clio/projects/proofs/code-1119-Q597/ct6.py').read().split('\ndef work(w):')[0])
def blockvec(X, Y):
    if len(X) != len(Y) or not X: return None
    v = tvec(X[0], Y[0])
    return v if v is not None and all(tr(s, v) == t for s, t in zip(X, Y)) else None
def cuts(mA, mB):
    out = []
    if len(mA) != len(mB): return ['len differ']
    for s in range(1, len(mA)):
        if blockvec(mA[s:], mB[:len(mA)-s]) is not None and blockvec(mA[:s], mB[len(mA)-s:]) is not None: out.append(s)
    return out
def work(w):
    st = Counter(); ns = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'C', DV)
                    b12, Jyx = mvl(b1, y, 'C', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                    X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                    K = set(X) | set(Y) | set(Yx) | set(Xy)
                    I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                    sides = {k[0] for k in I}
                    if len(sides) != 1: continue
                    side = sides.pop()
                    cl, ot = ((Jy, Jyx), (Jx, Jxy)) if side == 0 else ((Jx, Jxy), (Jy, Jyx))
                    for nm, (A, B) in (('clean', cl), ('other', ot)):
                        mA, mB, _, _ = middles(A, B)
                        c = cuts(mA, mB)
                        L = len(mA)
                        lab = 'len differ' if c == ['len differ'] else ('no cut' if not c else '+'.join('fwd(s=1)' if s == 1 else ('inv(s=L-1)' if s == L-1 else f's={s}') for s in c))
                        st[(side, nm, lab)] += 1
                        st[(side, nm, 'len', L, len(mB))] += 1
    return st, ns
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(NCH)))
    st = Counter(); n = 0
    for s, k in res: st.update(s); n += k
    print('SIZE', SIZE, 'states', n)
    for k in sorted(st, key=str): print(st[k], k)
