"""PROVE 2026-11-26 c2 B/B CONTROL (copy of code-1126-Q597/step0.py with B/B non-commuting pairs, x<y). Was: Step 0: shape of the other role E's divergence, per diverging B/C pair.  usage: python3 step0.py SIZE NW NCH [lo-hi]"""
import sys
CHARG = sys.argv[4] if len(sys.argv) > 4 else None
exec(open("/home/clio/projects/proofs/code-1125c2-Q597/bsall.py").read().split("\ndef work(w):")[0])
import pickle
def work(w):
    st = Counter(); recs = []
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    if not x < y: continue
                    b1, Jx = mvl(b, x, 'B', DV); b2, Jy = mvl(b, y, 'B', DV)
                    b12, Jyx = mvl(b1, y, 'B', DV); b21, Jxy = mvl(b2, x, 'B', DV)
                    if tiles_of(b12) == tiles_of(b21): continue
                    roles = [('B', Jx, Jxy), ('B', Jy, Jyx)]
                    info = []
                    for T, A, Bq in roles:
                        mA, mB, p, q = middles(A, Bq)
                        info.append(dict(T=T, A=A, Ar=Bq, p=p, tA=typed(A), tB=typed(Bq), mA=mA, mB=mB))
                    Ls = [r for r in (0, 1) if info[r]['p'] < len(info[r]['A']) and info[r]['tA'][info[r]['p']] == (1, aT[info[r]['T']])]
                    st[('#L', len(Ls))] += 1
                    if len(Ls) != 1: continue
                    L = Ls[0]; E = 1 - L
                    iL, iE = info[L], info[E]
                    pL, pE = iL['p'], iE['p']
                    case = 'b' if iL['tB'][pL] == (1, 'b') else 'a'
                    # shared square: tiles removed by both journeys at their contact steps
                    def removed(s):
                        return set(t for t in s[1] if t[0] == 'sq')
                    sqL = removed(iL['A'][pL]); sqE = removed(iE['A'][pE])
                    st[(case, 'E at p', ws([iE['tA'][pE]]), 'Er at p', ws([iE['tB'][pE]]), 'L at p', ws([iL['tA'][pL]]), 'Lr', ws([iL['tB'][pL]]))] += 1
                    st[(case, 'shared square at both contacts', bool(sqL & sqE))] += 1
                    st[(case, 'E mid len', len(iE['mA']), 'L mid len', len(iL['mA']))] += 1
                    st[(case, 'E midwords', ws(typed(iE['A'])[pE:pE+len(iE['mA'])]), '->', ws(typed(iE['Ar'])[pE:pE+len(iE['mB'])]))] += 1
                    recs.append((case, L, [ (i['T'], i['A'], i['Ar'], i['p'], len(i['mA']), len(i['mB'])) for i in info ]))
    pickle.dump(recs, open('/tmp/step0BB_%d_%d.pkl' % (SIZE, w), 'wb'))
    return st
if __name__ == '__main__':
    idx = list(range(*map(int, CHARG.split('-')))) if CHARG else list(range(NCH))
    tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for st in pool.imap_unordered(work, idx): tot.update(st)
    print('SIZE', SIZE, 'chunks', CHARG or 'all', 'of', NCH)
    for k in sorted(tot, key=str): print(tot[k], k)
