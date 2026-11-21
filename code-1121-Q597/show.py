"""PROVE 11-21: print a few clean->B and clean->C records: both journeys as (class,type, s0, s1) relative to Jo's start base,
up to and including the contact step, plus the contact tiles and which Jo step removed them.  usage: python3 show.py 5 NMAX"""
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
NMAX = int(os.environ.get('NMAX', 6))
def segs(s):
    vs = sorted(rh(s[1])); s0 = sorted({tuple(Fr(c) for c in v[:2]) for v in vs}); s1 = sorted({tuple(Fr(c) for c in v[2:]) for v in vs})
    return s0, s1
def rel(pts, o): return [tuple(int(a - b) for a, b in zip(P, o)) for P in pts]
shown = Counter()
for idx, (key, b) in enumerate(census_sized()):
    for b, DV in states(b, cap=CAP):
        M0 = bare(b); cs = corners(set(b.nest[0]['pos']))
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
                tag, A, B, Jo, cc, oc = (('clean->C', Jy, Jyx, Jx, y, x) if side == 0 else ('clean->B', Jx, Jxy, Jy, x, y))
                mA, mB, p, q = middles(A, B)
                if not bs_test(mA, mB): continue
                if shown[tag] >= NMAX: continue
                shown[tag] += 1
                o0, o1 = segs(Jo[0]); O = o0[0] + o1[0]
                print('=====', tag, 'clean corner', cc, 'Jo corner', oc, 'p', p, 'Jo base', O)
                Mx0 = apply(M0, Jo); W = set(M0) - set(Mx0)
                TA = typed(A); TO = typed(Jo)
                print(' clean:')
                for k, s in enumerate(A[:p+2]):
                    s0, s1 = segs(s); hit = [u[0] for u in s[1] if u in W]
                    print('   %2d %s s0 %s s1 %s %s' % (k, TA[k], rel(s0, O[:2]), rel(s1, O[2:]), ('CONTACT ' + str(hit)) if hit else ''))
                rem = {u: j for j, s in enumerate(Jo) for u in s[1]}
                print(' Jo:')
                for k, s in enumerate(Jo[:14]):
                    s0, s1 = segs(s)
                    hitk = [u[0] for u in A[p][1] if u in W and rem.get(u) == k]
                    print('   %2d %s s0 %s s1 %s %s' % (k, TO[k], rel(s0, O[:2]), rel(s1, O[2:]), ('<- removes contact ' + str(hitk)) if hitk else ''))
        if all(shown[t] >= NMAX for t in ('clean->B', 'clean->C')): raise SystemExit
