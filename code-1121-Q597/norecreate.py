"""PROVE 11-21: Lemma (no re-creation): in a terminating journey no tile removed by a step is created by a later step.
Part 1 (engine, dumps): direct check on the 11-20 journey dumps.
Part 2 (word model, every size): created-tile 4-offsets Q'_T, removed Q_T (as words.offs); window m = max |q - q'|_1;
then every word of length <= 2m+2 from each allowed class, both T: created_k  ∩  removed_j  = ∅ for j < k.
Control: created_k ∩ removed_j for j > k (a created tile removed LATER) -- excluded by Theorem main, must also be 0;
live control: created_k ∩ (removed_j translated by an edge vector of R_j) must fire.
usage: python3 norecreate.py NW"""
import sys, os, pickle, glob, itertools
from multiprocessing import Pool
from collections import Counter
sys.path.insert(0, '/home/clio/projects/proofs/code-1120-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597')
import words as Wd
def it(t): return (t[0], tuple(sorted(tuple(int(x) for x in P) for P in t[1])))
def tr(t, v): return (t[0], tuple(sorted(tuple(a+b for a, b in zip(P, v)) for P in t[1])))
TYPES = ('a0', 'a1', 'b')
def prep(T):
    C = Wd.configs(T); D = {}; Q, Qp = set(), set()
    for (ck, typ), (Rc, rem, Rn, cre) in C.items():
        for t in rem: Q.add(Wd.offs(t, Rc, T))
        for t in cre: Qp.add(Wd.offs(t, Rc, T))
        EV = Wd.edgev(Rc)
        ctl = frozenset(tr(it(t), tuple(int(x) for x in v)) for t in rem for v in EV if any(v))
        D[(ck, typ)] = (frozenset(it(t) for t in rem), frozenset(it(t) for t in cre), ctl, min(it(Rn)[1]), Wd.classkey(Rn))
    return D, Q, Qp
def job(args):
    T, ck0, pre, L = args
    D, _, _ = prep(T); st = Counter()
    def dfs(ck, b, hist, word):
        for typ in TYPES:
            if len(hist) < len(pre) and typ != pre[len(hist)]: continue
            rem, cre, ctl, bn, ckn = D[(ck, typ)]
            f = lambda S: frozenset(tr(t, b) for t in S)
            remk, crek, ctlk = f(rem), f(cre), f(ctl)
            for (rj, cj, tj) in hist:
                st['pairs'] += 1
                if crek & rj: st['RECREATION'] += 1
                if cj & remk: st['created then removed (Thm main)'] += 1
                if crek & tj: st['control fires'] += 1
            st['nodes'] += 1
            nh = hist + [(remk, crek, ctlk)]
            if len(nh) < L: dfs(ckn, tuple(x+y for x, y in zip(bn, b)), nh, word + [typ])
    dfs(ck0, (0, 0, 0, 0), [], [])
    return st
if __name__ == '__main__':
    NW = int(sys.argv[1])
    # Part 1: engine dumps
    src = open('/home/clio/projects/proofs/code-1118c2-Q597/bs6.py').read(); ns = {}
    exec(src[src.index('def pre(A, B):'):src.index('import migrate')], ns)
    e = Counter()
    for fn in glob.glob('/home/clio/projects/proofs/code-1120-Q597/dumpj5-*.pkl'):
        d = pickle.load(open(fn, 'rb'))
        for kind in ('BC', 'BB2'):
            for TT, M0, J in d[kind]:
                e['journeys'] += 1; removed = set()
                for s in J:
                    nr = {t for t in s[2] if t[0] != 'rh'}
                    if nr & removed: e['ENGINE RECREATION'] += 1
                    removed |= {t for t in s[1] if t[0] != 'rh'}
    print('Part 1 (engine dumps):', dict(e), flush=True)
    # Part 2: window and words
    jobs = []
    for T in 'BC':
        D, Q, Qp = prep(T)
        m = max(sum(abs(a - b) for a, b in zip(q, qp)) for q in Q for qp in Qp)
        L = 2*m + 2
        print('T', T, 'Q', sorted(Q), "Q'", sorted(Qp), 'm', m, 'window length L', L, flush=True)
        for ck in {k[0] for k in D}:
            for pre in itertools.product(TYPES, repeat=3): jobs.append((T, ck, pre, L))
    st = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for s in pool.imap_unordered(job, jobs): st.update(s)
    for k in sorted(st): print(st[k], k)
