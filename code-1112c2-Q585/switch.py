"""Split the computed identity  M∘β = β∘C  into two halves via BSS-style tableau switching Σ.

Σ(S): S an LR tableau on ν/κ (content ρ).  Put u_κ on κ.  Process the cells of κ in decreasing standard order of
u_κ (larger entry first; equal entries right-to-left).  For each such cell, slide S inward into it (ordinary jdt:
the hole takes the smaller of its right and next-row neighbours, ties -> next row).  The outer cell vacated by the
slide is labelled with the u_κ entry of the starting cell.  The labels form S* on ν/ρ, content κ.

Tests (all c>=2 cells as in compare.py, gl_4, |λ|,|μ|<=5, plus c>=3 up to 7):
  H1:  M(T) == Σ(T)            for T = β_CM(b)  (inputs on ν/μ)   -- Purbhoo half
  H2:  Σ(β_Q(b)) == β_Q(C(b))                                      -- commutor half, under the OTHER bijection
  H2': Σ(β_CM(b)) == β_CM(C(b))
  ctl: Σ is not lex/revlex (count cells where Σ differs from each), Σ∘Σ = id, and lesion Σ_bad
       (ties -> right neighbour, i.e. a broken slide) must fail.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1112-Q585')
from commutor import LR, C1, parts
from flock import is_lr, content
from compare import beta_cm, M, trim
from betaq import beta_q


def sigma(S, outer, tie='next', order='std'):
    S = dict(S); outer = list(outer)
    inner = [0]*len(outer)
    for (c, r) in S: pass
    for r, l in enumerate(outer):
        cs = [c for (c, rr) in S if rr == r]
        inner[r] = min(cs) if cs else l
    kappa = [x for x in inner if x]
    cells = [(c, r) for r, l in enumerate(kappa) for c in range(l)]
    if order == 'std':
        cells.sort(key=lambda x: (-(x[1]), -x[0]))          # u_κ entry = row+1 ; larger first, right-to-left
    else:
        cells.sort(key=lambda x: (-x[0], -(x[1])))          # lesion: column order
    out = {}
    for (c, r) in cells:
        hc, hr = c, r
        while True:
            R = S.get((hc+1, hr)); D = S.get((hc, hr+1))
            if R is None and D is None: break
            if D is not None and (R is None or D < R or (D == R and tie == 'next')):
                S[(hc, hr)] = D; del S[(hc, hr+1)]; hr += 1
            else:
                S[(hc, hr)] = R; del S[(hc+1, hr)]; hc += 1
        out[(hc, hr)] = r+1
    return out


def run(cells_spec):
    st = dict(cells=0, elems=0, H1=0, H2Q=0, H2CM=0, invol=0, lr_out=0, ctl_bad_differs=0,
              Sigma_ne_lex=0, stdcheck=0)
    for (lo_c, maxn) in cells_spec:
        for a in range(1, maxn+1):
            for bb in range(1, maxn+1):
                for lam in parts(a, 4):
                    for mu in parts(bb, 4):
                        L1 = LR(lam, mu, 4)
                        for nu, bs in L1.items():
                            if not (lo_c <= len(bs)) or (lo_c == 3 and a <= 5 and bb <= 5): continue
                            nt = trim(nu); st['cells'] += 1; bad_cell = False; lexd = False
                            L2 = sorted(LR(mu, lam, 4)[nu]); src = sorted(bs)
                            for b in bs:
                                st['elems'] += 1
                                cb = C1(b, lam, mu, 4)[0]
                                T = beta_cm(b, lam, mu, nt)
                                s = sigma(T, nt)
                                st['lr_out'] += is_lr(s) and content(s) == trim(mu)
                                st['invol'] += sigma(s, nt) == T
                                st['H1'] += M(T, nt) == s
                                st['H2CM'] += s == beta_cm(cb, mu, lam, nt)
                                st['H2Q'] += sigma(beta_q(b, lam, mu), nt) == beta_q(cb, mu, lam)
                                if sigma(T, nt, tie='right') != s or sigma(T, nt, order='col') != s: bad_cell = True
                                lexd |= s != beta_cm(L2[src.index(b)], mu, lam, nt)
                            st['ctl_bad_differs'] += bad_cell; st['Sigma_ne_lex'] += lexd
    return st


if __name__ == '__main__':
    t0 = time.time()
    print('c>=2, n<=5:', run([(2, 5)]))
    print('c>=3, 5<n<=7:', run([(3, 7)]), 'time %.0fs' % (time.time()-t0))
