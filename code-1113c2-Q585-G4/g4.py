"""G4: Purbhoo's Cor 5.3 map f -> g~ equals Sigma o K on Tao readouts, K = hw o revrot.

Phi (Purbhoo 0705.1184 p.20, steps 2-4, read with the typing fixed by the engine, see the .tex):
   f~  = rotate f by 180 deg, entries -> position in the REVERSE standard order  (a standard tableau)
   g   = straight LR tableau u_kappa of the inner shape kappa of f~
   g~  = slide the boxes of f~, in increasing order, through g (reverse jdt of g; ties -> the box BELOW)
Sigma_in (switch.sigma): slide f~ (or any semistandard X) INWARD into the cells of u_kappa, larger u-entry first.

Checks, per element of the census (engine: Tao(P), kappa(P)=T, Phi_G(P)):
  L0  revrot(tao) occupies exactly the cells of T (after translation)          [typing of Cor 5.3]
  L1  std(revrot(tao)) == f~ (rotate + reverse-standard positions)            [Lemma: std commutes with revrot]
  L2  vacated-cell sequence of Sigma_in identical for Y, std Y, hw Y          [shape evolution invariance]
  L3  coplactic: e_i commutes with every single slide along the Sigma_in run  (all i, all slides)
  L4  Phi_out(f~) == Sigma_in(f~)                                              [BSS order-independence, recalled]
  G4  Phi(tao) == Sigma(K(tao))       and  == Phi_G(P) (engine)  [the latter is H1 given K1]
  CTRL a: Sigma of the OTHER LR tableau(s) of the cell != Phi(tao)  (must fire on every element)
  CTRL b: Phi without the rotation (std(tao) in place) is not even the right shape / differs
  CTRL c: Phi_out with ties -> WEST (wrong) must break somewhere (else ties never met: report count)
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0, '/home/clio/projects/proofs/code-1113-Q585-H1')
sys.path.insert(0, '/home/clio/projects/proofs/code-1113c2-Q585-G4')
from collections import Counter
from explore import mosaic_of, elements
from readouts import readout, norm
from test_K import hw, revrot
from switch import sigma
from flock import is_lr, lr_tableaux
from falsifier import e_i

def std(X):
    """standardisation: order by (value, column) -- equal entries west first (Purbhoo's standard order)."""
    order = sorted(X, key=lambda x: (X[x], x[0]))
    return {x: k + 1 for k, x in enumerate(order)}

def rot_revstd(f, Cc, Rr):
    """Purbhoo step 2 on a tableau placed in a box with columns 0..Cc-1, rows 0..Rr-1."""
    order = sorted(f, key=lambda x: (f[x], x[0]), reverse=True)  # reverse standard order
    return {(Cc - 1 - c, Rr - 1 - r): k + 1 for k, (c, r) in enumerate(order)}

def inner_of(X, outer):
    """inner shape of X inside the straight shape outer (rows with no X-cell are entirely inner)."""
    return [min([c for (c, rr) in X if rr == r], default=l) for r, l in enumerate(outer)]

def sigma_trace(X, outer):
    """Sigma_in with the vacated-cell sequence, and the full list of intermediate X's (for L3)."""
    X = dict(X); kap = [x for x in inner_of(X, outer) if x]
    cells = sorted(((c, r) for r, l in enumerate(kap) for c in range(l)), key=lambda x: (-x[1], -x[0]))
    vac = []; states = [dict(X)]
    for (c, r) in cells:
        hc, hr = c, r
        while True:
            R = X.get((hc + 1, hr)); D = X.get((hc, hr + 1))
            if R is None and D is None: break
            if D is not None and (R is None or D <= R):
                X[(hc, hr)] = D; del X[(hc, hr + 1)]; hr += 1
            else:
                X[(hc, hr)] = R; del X[(hc + 1, hr)]; hc += 1
        vac.append((hc, hr)); states.append(dict(X))
    out = {v: cells[k][1] + 1 for k, v in enumerate(vac)}
    return out, vac, states, cells

def one_slide(X, cell):
    X = dict(X); hc, hr = cell
    while True:
        R = X.get((hc + 1, hr)); D = X.get((hc, hr + 1))
        if R is None and D is None: break
        if D is not None and (R is None or D <= R):
            X[(hc, hr)] = D; del X[(hc, hr + 1)]; hr += 1
        else:
            X[(hc, hr)] = R; del X[(hc + 1, hr)]; hc += 1
    return X

TIES = Counter()
def phi_out(ft, outer, tie='south'):
    """Purbhoo step 4: g = u_kappa on the inner shape; for k = 1..N slide the hole at f~^{-1}(k) back through g."""
    kap = [x for x in inner_of(ft, outer) if x]
    g = {(c, r): r + 1 for r, l in enumerate(kap) for c in range(l)}
    for k in range(1, len(ft) + 1):
        x = next(y for y, v in ft.items() if v == k)
        hc, hr = x
        while True:
            W = g.get((hc - 1, hr)); S = g.get((hc, hr - 1))
            if W is None and S is None: break
            if W is not None and S is not None and W == S: TIES['met'] += 1
            takeS = S is not None and (W is None or S > W or (S == W and tie == 'south'))
            if takeS:
                g[(hc, hr)] = S; del g[(hc, hr - 1)]; hr -= 1
            else:
                g[(hc, hr)] = W; del g[(hc - 1, hr)]; hc -= 1
        # f~ entry k now sits at (hc,hr): irrelevant for g~
    return g

def Y_unrot(tao, Y):
    """CTRL b: tao translated onto Y's cells WITHOUT rotating (only meaningful if shapes coincide; else tao as is)"""
    Z = translate_onto(tao, set(Y)); return Z if Z is not None else Y_bad(Y)
def Y_bad(Y):
    # shape does not even match: control fires trivially; return the row-reading standardisation reversed
    order = sorted(Y, key=lambda x: (-x[1], x[0])); return {x: k + 1 for k, x in enumerate(order)}

def translate_onto(Y, target_cells):
    """translate Y so that its cells coincide with target_cells (both normalised); None if impossible"""
    ny = dict(norm(Y)); c0 = min(c for c, r in target_cells); r0 = min(r for c, r in target_cells)
    Z = {(c + c0, r + r0): v for (c, r), v in ny.items()}
    return Z if set(Z) == set(target_cells) else None

def run(amax, bmax, cmin, label):
    st = Counter(); TIES.clear()
    for lam, mu, nt, T in elements(amax, bmax, cmin):
        st['n'] += 1; l = len([x for x in lam if x])
        b, k0, k1 = mosaic_of(T, nt)
        _, tao = readout(b, 0, 1, 'N'); _, PG = readout(b, 1, -1, 'W')
        tao = dict(norm(tao))
        Y = translate_onto(revrot(tao, l), set(T)); st['L0'] += Y is not None
        if Y is None: continue
        # place tao in a box so that its 180-rotation lands on Y's cells: box = bounding box of tao + Y offsets
        Cc = max(c for c, r in Y) + min(c for c, r in tao) + 1 + 0
        # direct: rotation (c,r)->(Cc-1-c, Rr-1-r) must send tao's cells onto Y's cells
        cY = max(c for c, r in Y) + min(c for c, r in tao); rY = max(r for c, r in Y) + min(r for c, r in tao)
        ft = rot_revstd(tao, cY + 1, rY + 1)
        st['L1 cells'] += set(ft) == set(Y)
        st['L1'] += ft == std(Y)
        H = hw(Y)
        s1, v1, S1, cells = sigma_trace(Y, nt); s2, v2, _, _ = sigma_trace(std(Y), nt); s3, v3, _, _ = sigma_trace(H, nt)
        st['L2'] += (v1 == v2 == v3)
        # L3 coplactic along the run on Y: for each intermediate state and each i, e_i(slide(Z)) == slide(e_i(Z))
        ok3 = True
        for k, cell in enumerate(cells):
            Z = S1[k]
            for i in range(1, l):
                ez = e_i(Z, i)
                if ez is None: continue
                st['L3 tests'] += 1
                if e_i(one_slide(Z, cell), i) != one_slide(ez, cell): ok3 = False
        st['L3'] += ok3
        po = phi_out(ft, nt); st['L4'] += norm(po) == norm(sigma(ft, nt))
        g4 = norm(po) == norm(sigma(H, nt)); st['G4 Phi=Sigma.K'] += g4
        st['Phi = Phi_G (engine)'] += norm(po) == norm(PG)
        st['Sigma.K is LR'] += is_lr(sigma(H, nt))
        others = [O for O in lr_tableaux(set(T), tuple(x for x in lam if x)) if O != H]
        st['CTRL a fires'] += bool(others) and all(norm(sigma(O, nt)) != norm(po) for O in others)
        st['CTRL b shape-matched'] += translate_onto(tao, set(Y)) is not None
        st['CTRL b fires'] += norm(phi_out(std(tao), nt)) if False else norm(phi_out(std(Y_unrot(tao, Y)), nt)) != norm(po)
        st['CTRL c fires'] += norm(phi_out(ft, nt, tie='west')) != norm(po)
    print('==', label, dict(st), 'ties met in phi_out:', TIES['met'])

if __name__ == '__main__':
    run(4, 4, 2, 'c>=2, |lam|,|mu|<=4')
    run(5, 5, 2, 'c>=2, |lam|,|mu|<=5')
