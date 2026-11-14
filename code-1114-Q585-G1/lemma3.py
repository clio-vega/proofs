"""Checks for the reduction lemmas.
 (L-flock) unique LR flock on alpha in + frame = row index; on the 180-rotated alpha (- frame) = alpha'_c - r ('depth').
           Pure check via lr_tableaux for all alpha, |alpha|<=9 (independent of the migration engine).
 (L3) on the c>=2 census: arrival sequence of (ii) = DEcreasing standard order of R_ii; of (iii) = INcreasing std order of R_iii;
      reversal c_k = b_{m+1-k} (after normalisation).  Control: unreversed c_k = b_k must FAIL somewhere."""
import sys, copy
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0,'/home/clio/projects/proofs/code-1113-Q585-H1')
from readouts import *
from compare import beta_cm, trim
from commutor import LR, parts
from collections import Counter
if __name__ == '__main__':
    st = Counter()
    for size in range(1, 10):
        for al in parts(size, size):
            al = tuple(x for x in al if x)
            cells = {(c, r) for r, l in enumerate(al) for c in range(l)}
            conj = [sum(1 for l in al if l > c) for c in range(al[0])]
            L = lr_tableaux(cells, al); st['+ unique'] += len(L) == 1
            st['+ = row'] += L[0] == {(c, r): r+1 for c, r in cells}
            R0 = len(al) - 1; C0 = al[0] - 1
            rot = {(C0-c, R0-r) for c, r in cells}
            L2 = lr_tableaux(rot, al); st['- unique'] += len(L2) == 1
            st['- = depth'] += L2[0] == {(C0-c, R0-r): conj[c]-r for c, r in cells}
            st['alphas'] += 1
    print('flock lemma', dict(st)); st = Counter()

def std_key(T):  # standard order: smaller value, or equal and west
    return sorted(T, key=lambda b: (T[b], b[0]))
def arrivals(b, sgn, dirn):
    u, w = b.nest[0]['u'], b.nest[0]['w']; E_ = tuple(sgn*a for a in u); N_ = tuple(sgn*a for a in w)
    fl = unique_flock(b, 0, E_, N_)
    dvec = {'E': E_, 'N': N_, 'W': tuple(-a for a in E_), 'S': tuple(-a for a in N_)}[dirn]
    key = {t: (v, box_of(b, 0, b.tiles[t][1], E_, N_)[0]) for t, v in fl.items()}
    order = sorted(fl, key=lambda t: key[t], reverse=(dirn in 'NE'))
    inv = {t: ij for ij, t in b.nest[0]['pos'].items()}; land = []
    for t in order:
        del b.nest[0]['pos'][inv[t]]; R, X, ij, s = b.journey(t, 0, unit_dir(dvec)); land.append((R, fl[t]))
    E2, N2, _ = induced(b, X, E_, N_)
    boxes = [box_of(b, X, b.tiles[R][1], E2, N2) for R, v in land]
    c0 = min(c for c, r in boxes); r0 = min(r for c, r in boxes)
    boxes = [(c-c0, r-r0) for c, r in boxes]
    return boxes, [v for R, v in land]
if __name__ == '__main__':
    for a in range(1, 6):
     for c_ in range(1, 6):
      for lam in parts(a, 4):
       for mu in parts(c_, 4):
        for nu, bs in LR(lam, mu, 4).items():
          if len(bs) < 2: continue
          nt = trim(nu)
          for bb in bs:
            T = beta_cm(bb, lam, mu, nt)
            b0 = canonical(nt, max(nt[0], 1)+1, len(nt)+1); B = b0.nest[1]
            migrate_flock(b0, 1, {B['pos'][(c, r)]: v for (c, r), v in T.items()}, B['u'], B['w'], 'E')
            bi, vi = arrivals(copy.deepcopy(b0), 1, 'N'); bj, vj = arrivals(copy.deepcopy(b0), -1, 'W')
            Ti = dict(zip(bi, vi)); Tj = dict(zip(bj, vj)); st['n'] += 1
            st['(ii) arrivals = decreasing std'] += bi == std_key(Ti)[::-1]
            st['(iii) arrivals = increasing std'] += bj == std_key(Tj)
            st['values (ii) reversed = values (iii)'] += vi[::-1] == vj
            st['REVERSAL c_k = b_{m+1-k}'] += bj == bi[::-1]
            st['control c_k = b_k holds'] += bj == bi
            st['R_ii = R_iii'] += Ti == Tj
    print(dict(st))
