"""G4 as a pure tableau identity, no mosaics: for every LR tableau f on A/B inside an R x C box (French, (c,r)),
    Phi(f) := Purbhoo outward switching of u_{A^vee} through f~ = rot180(reverse-standard positions of f)
    equals Sigma(revrot f) and Sigma(hw revrot f)   (Sigma = inward switching, switch.sigma).
Also: Sigma(X) is constant on the whole crystal component of revrot f (all elements, not just hw) -- the
coplactic invariance used in the proof -- and Phi is a bijection LR(A/B;c) -> LR(B^vee/c ; A^vee).
Control: Sigma of a DIFFERENT LR tableau of the same type must differ (c>=2 types)."""
import sys, itertools, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585'); sys.path.insert(0, '/home/clio/projects/proofs/code-1113-Q585-H1')
sys.path.insert(0, '/home/clio/projects/proofs/code-1113c2-Q585-G4')
from collections import Counter
from flock import lr_tableaux, is_lr, content
from switch import sigma
from test_K import hw, revrot
from g4 import phi_out, rot_revstd, std
from falsifier import e_i, f_i

def partitions_in_box(R, C):
    for t in itertools.product(range(C + 1), repeat=R):
        if all(t[i] >= t[i + 1] for i in range(R - 1)): yield t

def parts_of(n, k):
    def rec(n, mx, k):
        if n == 0: yield (); return
        if k == 0: return
        for x in range(min(n, mx), 0, -1):
            for rest in rec(n - x, x, k - 1): yield (x,) + rest
    yield from rec(n, n, k)

def component(X, n):
    seen = {tuple(sorted(X.items()))}; todo = [X]
    while todo:
        Y = todo.pop()
        for i in range(1, n):
            for op in (e_i, f_i):
                Z = op(Y, i)
                if Z is not None:
                    k = tuple(sorted(Z.items()))
                    if k not in seen: seen.add(k); todo.append(Z)
    return [dict(k) for k in seen]

def run(R, C, nmax, gl):
    st = Counter(); t0 = time.time()
    for A in partitions_in_box(R, C):
        for B in partitions_in_box(R, C):
            if any(B[i] > A[i] for i in range(R)): continue
            cells = {(c, r) for r in range(R) for c in range(B[r], A[r])}
            if not cells or len(cells) > nmax: continue
            Bv = tuple(C - B[R - 1 - r] for r in range(R))
            for cont in parts_of(len(cells), gl):
                Ls = lr_tableaux(cells, cont)
                if not Ls: continue
                st['types'] += 1; st['c>=2 types'] += len(Ls) >= 2
                outs = []
                for f in Ls:
                    st['f'] += 1; l = len(cont)
                    Y = {(C - 1 - c, R - 1 - r): l + 1 - v for (c, r), v in f.items()}   # revrot inside the box
                    ft = rot_revstd(f, C, R)
                    st['std revrot = f~'] += std(Y) == ft
                    p = phi_out(ft, Bv); s0 = sigma(Y, Bv); s1 = sigma(hw(Y), Bv)
                    st['Phi = Sigma(revrot f)'] += p == s0
                    st['Phi = Sigma(hw revrot f)'] += p == s1
                    st['Phi is LR'] += is_lr(p)
                    comp = component(Y, max(l, 2)); st['component elements'] += len(comp)
                    st['Sigma constant on component'] += all(sigma(Z, Bv) == s0 for Z in comp)
                    outs.append(tuple(sorted(p.items())))
                st['Phi injective on type'] += len(set(outs)) == len(outs)
                if len(Ls) >= 2:
                    st['CTRL distinct inputs give distinct Sigma'] += len({tuple(sorted(sigma(hw({(C - 1 - c, R - 1 - r): len(cont) + 1 - v for (c, r), v in f.items()}), Bv).items())) for f in Ls}) == len(Ls)
    print(f'box {R}x{C}, |A/B|<={nmax}, gl_{gl}:', dict(st), '%.0fs' % (time.time() - t0))

if __name__ == '__main__':
    run(3, 3, 9, 3)
    run(4, 4, 9, 4)
    run(4, 5, 10, 4)
