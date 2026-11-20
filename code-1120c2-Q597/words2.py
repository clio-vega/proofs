"""PROVE 11-20 c2: Lemma AW, fast exhaustive version of words.py (integer points, set intersections, Pool over prefixes).
For every word of length <= L in {a0,a1,b} from each allowed class, both T, and every pair of steps k != j in the word:
  AW:      (ascending R_j-edge translates of tiles removed at j)  ∩  (tiles removed at k)  = ∅.
  CONTROL: (ascending NON-edge unit translates of tiles removed at j) ∩ (tiles removed at k)  -- must be nonempty sometimes.
Words from words.configs (asserted unique per (class,type), sign-consistent).   usage: python3 words2.py L NW"""
import sys, itertools
from multiprocessing import Pool
from collections import Counter
import words as Wd
from psilib import psi
def it(t): return (t[0], tuple(sorted(tuple(int(x) for x in P) for P in t[1])))
def tr(t, v): return (t[0], tuple(sorted(tuple(a+b for a, b in zip(P, v)) for P in t[1])))
def prep(T):
    C = Wd.configs(T); D = {}
    for (ck, typ), (Rc, rem, Rn, cre) in C.items():
        R = it(Rc); assert min(R[1]) == (0, 0, 0, 0)
        EV = Wd.edgev(Rc)
        asc, ctl = set(), set()
        for u in rem:
            pu = psi(u, T); ui = it(u)
            for v in Wd.W12:
                c = 1 if any(v[:2]) else 0
                if psi(Wd.trs(u, v), T)[c] - pu[c] != 1: continue
                (asc if v in EV else ctl).add(tr(ui, v))
        Rn_i = it(Rn)
        D[(ck, typ)] = (frozenset(it(t) for t in rem), frozenset(asc), frozenset(ctl), min(Rn_i[1]), Wd.classkey(Rn))
    return D
TYPES = ('a0', 'a1', 'b')
def job(args):
    T, ck0, prefix, L = args
    D = prep(T); st = Counter()
    def step(ck, b, typ):
        rem, asc, ctl, bn, ckn = D[(ck, typ)]
        sh = b
        f = lambda S: frozenset(tr(t, sh) for t in S)
        return f(rem), f(asc), f(ctl), tuple(x+y for x, y in zip(bn, sh)), ckn
    # replay prefix
    hist = []; ck, b = ck0, (0, 0, 0, 0)
    def dfs(ck, b, hist, depth_word):
        for typ in TYPES:
            if depth_word is not None and len(hist) < len(depth_word) and typ != depth_word[len(hist)]: continue
            rem, asc, ctl, bn, ckn = step(ck, b, typ)
            for (rk, ak, ck_) in hist:
                if asc & rk or ak & rem: st['AW VIOLATION'] += 1
                if ctl & rk or ck_ & rem: st['control fires'] += 1
                st['pairs checked'] += 1
            st['nodes'] += 1
            nh = hist + [(rem, asc, ctl)]
            if len(nh) < L: dfs(ckn, bn, nh, depth_word)
    dfs(ck0, (0, 0, 0, 0), [], prefix)
    return st
if __name__ == '__main__':
    L = int(sys.argv[1]); NW = int(sys.argv[2])
    jobs = []
    for T in 'BC':
        D = prep(T)
        for ck in {k[0] for k in D}:
            for pre in itertools.product(TYPES, repeat=3):
                jobs.append((T, ck, pre, L))
    st = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for s in pool.imap_unordered(job, jobs): st.update(s)
    print('L', L, 'jobs', len(jobs))
    for k in sorted(st, key=str): print(st[k], k)
