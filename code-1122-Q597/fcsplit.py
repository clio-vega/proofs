"""PROVE 11-22 Step 0: split (FC<=) -- does the divergence step read a W tile, and WHY.
Per record (clean B/C, 521 at <=5) and per role (BBy, BBx controls):
  A  = clean journey on M0 (toward T),  B = its re-run on Mx = M0.Jo,  Jo = the other journey on M0 (toward T').
  X_k = M0.A[:k],  Y_k = Mx.B[:k]  (k <= p, the common prefix).  W = M0 - Mx,  Xadd = Mx - M0.
Checks
  PRE   records 521; p = longest common prefix; fc (first A step removing a W tile) == p.
  DINV  X_k - Y_k == W and Y_k - X_k == Xadd for every k <= p   (Lemma Dinv, by hand)
  UNIQ  at every step of A, B and Jo: the good configurations of the target at R_k whose 3 non-rh tiles are all on the
        board are exactly one, and it is the step taken   (Lemma U, from Prop good / Lemma table of the no-repeat tex)
  P     at p: H_p ∩ W != 0;  H^B_p ∩ Xadd != 0;  H_p ⊄ Y_p;  H^B_p ⊄ X_p.
  WIN   (brief's table) the first k whose read window (the 9 non-rh tiles of the 3 good configs at R_k) meets
        W-on-X_k or Xadd-on-Y_k: kind (W only / Xadd only / both) x (k == p ?).
  CTRL  steps k < p whose window meets W ∪ Xadd while the step agrees (window difference without divergence), and steps
        where the SET of present good configs differs between X_k and Y_k (must be 0 before p by Dinv+U).
usage: KEEP=k python3 fcsplit.py SIZE NW NCH"""
import os, sys
exec(open("/home/clio/projects/proofs/code-1121-Q597/diverge.py").read().split("\ndef work(w):")[0])
_here = os.getcwd()
os.chdir('/home/clio/projects/proofs/code-1120c2-Q597'); sys.path.insert(0, '/home/clio/projects/proofs/code-1120c2-Q597')
exec(open('/home/clio/projects/proofs/code-1120c2-Q597/words.py').read().split('\ndef run(T, L, st):')[0])
os.chdir(_here)
CF = {T: configs(T) for T in 'BC'}
TYPES = ('a0', 'a1', 'b')
def goodcfgs(T, R):
    """type -> frozenset of the 3 non-rh tiles of the good T-configuration of that type at rhombus R (engine tile)."""
    Rn = norm(R); ck = classkey(Rn); out = {}
    for typ in TYPES:
        Rc, rem, Rnew, cre = CF[T][(ck, typ)]
        sh = base(Rn); assert trs(Rc, sh) == Rn
        out[typ] = frozenset(trs(t, sh) for t in rem)
    return out
def nn(S): return {norm(t) for t in S}
def present(T, R, board):
    return {typ for typ, G in goodcfgs(T, R).items() if G <= board}
def stepcheck(T, J, M, st, tag, who):
    Bd = nn(M)
    for k, s in enumerate(J):
        R = norm(rhof(s)); pres = present(T, R, Bd)
        g = {norm(t) for t in s[1] if t[0] != 'rh'}
        taken = [typ for typ, G in goodcfgs(T, R).items() if G == g]
        st[(tag, 'UNIQ', who, 'n present', len(pres), 'taken is the present one', len(taken) == 1 and pres == set(taken))] += 1
        Bd -= nn(s[1]); Bd |= nn(s[2])
def work(w):
    st = Counter(); nstates = 0
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            nstates += 1; M0 = bare(b)
            cs = corners(set(b.nest[0]['pos']))
            for x in cs:
                for y in cs:
                    if x == y: continue
                    for kind in ('BC', 'BB'):
                        if kind == 'BB' and not x < y: continue
                        tx, ty = kind
                        b1, Jx = mvl(b, x, tx, DV); b2, Jy = mvl(b, y, ty, DV)
                        b12, Jyx = mvl(b1, y, ty, DV); b21, Jxy = mvl(b2, x, tx, DV)
                        if kind == 'BB':
                            if tiles_of(b12) == tiles_of(b21): continue
                            roles = (('BBy', Jy, Jyx, Jx, 'B', 'B'), ('BBx', Jx, Jxy, Jy, 'B', 'B'))
                        else:
                            P, P1, P2, P12, P21 = map(pq, (b, b1, b2, b12, b21))
                            X, Y, Yx, Xy = dl(P, P1), dl(P, P2), dl(P1, P12), dl(P2, P21)
                            K = set(X) | set(Y) | set(Yx) | set(Xy)
                            I = {k for k in K if not (Yx.get(k, Z) == Y.get(k, Z) and Xy.get(k, Z) == X.get(k, Z))}
                            sides = {k[0] for k in I}
                            if len(sides) != 1: continue
                            side = sides.pop()
                            roles = (('clean->C', Jy, Jyx, Jx, 'C', 'B'),) if side == 0 else (('clean->B', Jx, Jxy, Jy, 'B', 'C'),)
                        for tag, A, B, Jo, T, Tp in roles:
                            mA, mB, p, q = middles(A, B)
                            if kind == 'BC' and not bs_test(mA, mB): continue
                            if len(mA) < 2 or len(mB) < 1: continue
                            st[(tag, 'records')] += 1
                            M0n = nn(M0); Mx = nn(apply(M0, Jo))
                            W = M0n - Mx; Xadd = Mx - M0n
                            # PRE: p is the longest common prefix and both journeys step at p
                            st[(tag, 'PRE p = lcp, A[p] != B[p]', all(A[k] == B[k] for k in range(p)) and A[p] != B[p])] += 1
                            rd = [k for k, s in enumerate(A) if nn(s[1]) & W]
                            st[(tag, 'PRE fc == p', bool(rd) and rd[0] == p)] += 1
                            # UNIQ on A (on M0), B (on Mx), Jo (on M0)
                            stepcheck(T, A, M0, st, tag, 'A'); stepcheck(T, B, apply(M0, Jo), st, tag, 'B(re-run)'); stepcheck(Tp, Jo, M0, st, tag, 'Jo')
                            # DINV along the prefix, WIN, CTRL
                            Xk, Yk = set(M0n), set(Mx); first = None; dinv = True
                            for k in range(p + 1):
                                dinv &= (Xk - Yk == W) and (Yk - Xk == Xadd)
                                R = norm(rhof(A[k]))
                                win = set().union(*goodcfgs(T, R).values())
                                mw = bool(win & W & Xk); mx = bool(win & Xadd & Yk)
                                if (mw or mx) and first is None:
                                    first = (k, mw, mx)
                                if k < p:
                                    st[(tag, 'CTRL k<p: window meets W/Xadd', (mw, mx), 'present sets differ', present(T, R, Xk) != present(T, R, Yk))] += 1
                                    Xk -= nn(A[k][1]); Xk |= nn(A[k][2]); Yk -= nn(B[k][1]); Yk |= nn(B[k][2])
                            st[(tag, 'DINV X_k-Y_k == W, Y_k-X_k == Xadd for all k<=p', dinv)] += 1
                            fk, mw, mx = first
                            st[(tag, 'WIN first window contact: W', mw, 'Xadd', mx, 'at k == p', fk == p)] += 1
                            Hp = nn(A[p][1]) - {norm(rhof(A[p]))}; HBp = nn(B[p][1]) - {norm(rhof(B[p]))}
                            st[(tag, 'P: Hp&W', len(Hp & W), 'HBp&Xadd', len(HBp & Xadd), 'Hp<=Y_p', Hp <= Yk, 'HBp<=X_p', HBp <= Xk)] += 1
                            st[(tag, 'LESION present on X_p | Y_p (union board; must be 2)', len(present(T, norm(rhof(A[p])), Xk | Yk)))] += 1
                            st[(tag, 'P: present on X_p', tuple(sorted(present(T, norm(rhof(A[p])), Xk))), 'on Y_p', tuple(sorted(present(T, norm(rhof(A[p])), Yk))))] += 1
    return st, nstates
if __name__ == '__main__':
    with Pool(NW, maxtasksperchild=1) as pool: res = list(pool.imap_unordered(work, range(int(os.environ.get('KEEP', NCH)))))
    st = Counter(); ns = 0
    for s, n in res: st.update(s); ns += n
    print('SIZE', SIZE, 'states', ns, 'KEEP', os.environ.get('KEEP', NCH), 'of', NCH)
    for k in sorted(st, key=str): print(st[k], k)
