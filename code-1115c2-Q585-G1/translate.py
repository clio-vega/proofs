"""Per-step test of the translation conjecture: source of (iii) step k on M == phi(source of (iii) step k on M_1), where
phi shifts cells of alpha-hat lying under the top row up by one row. Row coordinate determined from the top-row cells."""
import pickle
from collections import Counter
R = pickle.load(open('/tmp/c2-exchange.pkl','rb')); st = Counter(); bad = None
for r in R:
    top = r['top']
    if len(set(j for i, j in top)) == 1: row = 1   # (i,j): j is the row index
    elif len(set(i for i, j in top)) == 1: row = 0
    else: raise SystemExit('top row not a row')
    if len(top) == 1: st['top row has 1 cell: row coord ambiguous'] += 1
    rowof = lambda p: p[row]; colof = lambda p: p[1-row]
    toprow = rowof(top[0]); topcols = {colof(p) for p in top}
    def phi(p):
        if colof(p) in topcols: return (p[0], p[1]+1) if row == 1 else (p[0]+1, p[1])
        return p
    srcM = [s for s, l, v in r['M_iii']]; srcM1 = [s for s, l, v in r['M1_iii']]
    k = len(srcM1); st['n'] += 1
    t_ok = all(srcM[i] == phi(srcM1[i]) for i in range(k)); st['sources: s_k = phi(s1_k) all k'] += t_ok
    st['steps where phi is not the identity (non-trivial translation)'] += sum(phi(s) != s for s in srcM1)
    st['steps total'] += k
    u_ok = all(srcM[i] == srcM1[i] for i in range(k)); st['control: s_k = s1_k all k'] += u_ok
    # values: (iii) flock values (depths) agree?
    st['flock values agree'] += [v for s, l, v in r['M_iii']][:k] == [v for s, l, v in r['M1_iii']]
    st['last alpha_top steps of M move the bottom cells of top columns'] += {colof(s) for s in srcM[k:]} == topcols and all(rowof(s) == 0 for s in srcM[k:])
    if not t_ok and bad is None: bad = (r['key'], r['A'], top, srcM, srcM1)
print(dict(st)); print(bad)
