"""PROVE 11-20 c2: the abstract-word model of a journey, and Lemma AW (G1c at every size, finite check).
Fact (Lemma table + Prop good, 11-20): every step of a terminating journey is the UNIQUE good configuration of its
(class of R, type) -- checked here (assert).  So a journey is a word in {a0, a1, b}, and R_k and the tiles removed at
step k are determined by the word and R_0, up to translation.
Checks (board-free, exhaustive over all words of length <= L from each allowed class, both T):
  (AW)  if u is removed at step j and v is an edge vector of R_j with psi(u+v) = psi(u) + e_c (ascending), then u+v is
        removed at NO step k != j of the word.
  (Q)   offsets om(t,R) of removed tiles lie in Q_T; (LAM) Delta lambda in {0,1}^2, Delta mu in {0,-1}^2;
  (ST)  the steps with Delta(mu,lambda) = 0 are exactly V-type b-steps, and no two such steps are consecutive.
Control (must FIRE): the same test with v a NON-edge unit vector that is ascending.
usage: python3 words.py L"""
import sys
sys.path.insert(0, '../code-1120-Q597')
from fractions import Fraction as Fr
from collections import Counter
import localenum as LE
from psilib import psi, ETA
UN = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
W12 = [u + (0, 0) for u in UN] + [(0, 0) + u for u in UN]
def trs(t, v): return (t[0], frozenset(tuple(Fr(a)+b for a, b in zip(P, v)) for P in t[1]))
def norm(t): return (t[0], frozenset(tuple(Fr(x) for x in P) for P in t[1]))
def base(R): return min(R[1])
def configs(T):
    """(class key, type) -> (R canonical, removed tiles, R' ), all relative to base(R) = 0."""
    C = {}
    for u in LE.UNITS:
        if LE.ETA[T](u) == 0: continue
        R, Hs = LE.classify(T, u)
        for typ, H in Hs:
            r = LE.analyse(T, R, typ, set(H))
            if not r['xrule'] or r['k'] is None or r['across'] is False: continue
            V = set().union(*[t[1] for t in H]); zs = LE.convex_hull([LE.phi(P) for P in V])
            corners = [P for P in V if LE.is_corner(LE.phi(P), zs)]
            c = tuple(sum(Fr(P[i]) for P in corners)/6 for i in range(4))
            Hn = frozenset((t[0], frozenset(LE.rot_about(P, c, r['k']) for P in t[1])) for t in H)
            Rn = next(t for t in Hn if t[0] == 'rh'); R0 = norm(R)
            b0 = base(R0); sh = tuple(-x for x in b0)
            key = (frozenset({u, tuple(-x for x in u)}), typ)
            val = (trs(R0, sh), frozenset(trs(norm(t), sh) for t in H if t[0] != 'rh'), trs(norm(Rn), sh), frozenset(trs(norm(t), sh) for t in Hn if t[0] != 'rh'))
            if key in C: assert C[key] == val, ('sign-inconsistent config', T, key)
            C[key] = val
    return C
def classkey(R):
    P = sorted(R[1]); d = set()
    for A in P:
        for B in P:
            w = tuple(b-a for a, b in zip(A, B))
            if any(w[:2]) and not any(w[2:]): d.add(tuple(int(x) for x in w[:2]))
    assert len(d) == 2, d
    return frozenset(d)
def offs(t, R, T):
    pt, pr = psi(t, T), psi(R, T)
    at = (min(P[0] for P in t[1]), min(P[2] for P in t[1])); ar = (min(P[0] for P in R[1]), min(P[2] for P in R[1]))
    return (int(pt[0]-pr[0]), int(pt[1]-pr[1]), int(at[0]-ar[0]), int(at[1]-ar[1]))
def edgev(R):
    return {tuple(b-a for a, b in zip(A, B)) for A in R[1] for B in R[1]} - {(0, 0, 0, 0)}
def run(T, L, st):
    C = configs(T); types = sorted({k[1] for k in C})
    st[(T, 'configs', len(C))] += 1
    for key in {k[0] for k in C}:
        R0 = C[(key, types[0])][0]
        # DFS over words; stack entries: (R, history list of (R_k, removed_k))
        def dfs(R, hist, prevstat):
            if len(hist) == L: return
            ck = classkey(R)
            for typ in types:
                if (ck, typ) not in C: st[(T, 'missing config', typ)] += 1; continue
                Rc, rem, Rn, cre = C[(ck, typ)]
                sh = base(R); assert trs(Rc, sh) == R
                remk = frozenset(trs(t, sh) for t in rem); Rnew = trs(Rn, sh)
                j = len(hist)
                # (Q), (LAM), (ST)
                for t in remk: st[(T, 'Q offset', offs(t, R, T))] += 1
                d = tuple(a - b for a, b in zip(offs(Rnew, R, T), (0, 0, 0, 0)))
                st[(T, 'delta (mu0,mu1,lam0,lam1)', typ, d)] += 1
                stat = d == (0, 0, 0, 0)
                if stat and prevstat: st[(T, 'TWO CONSECUTIVE STATIONARY STEPS')] += 1
                # (AW): new step j vs all previous k, both roles
                EVj = edgev(R)
                newhist = hist + [(R, remk)]
                for k, (Rk, remkk) in enumerate(hist):
                    for (Ru, U, X, lab) in ((R, remk, remkk, 'u at later step j, u+v at earlier k'), (Rk, remkk, remk, 'u at earlier k, u+v at later j')):
                        EV = edgev(Ru)
                        for u in U:
                            pu = psi(u, T)
                            for v in W12:
                                c = 1 if any(v[:2]) else 0
                                w = trs(u, v)
                                if psi(w, T)[c] - pu[c] != 1: continue
                                isedge = v in EV
                                st[(T, 'AW' if isedge else 'CONTROL non-edge ascending', lab, 'u+v removed', w in X)] += 1
                dfs(Rnew, newhist, stat)
        dfs(R0, [], False)
if __name__ == '__main__':
    L = int(sys.argv[1]); st = Counter()
    for T in 'BC': run(T, L, st)
    for k in sorted(st, key=str): print(st[k], k)
