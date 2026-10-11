"""Q417 STEP 3: the MECHANISM, tested in general -- plus the one-line screen that
the mechanism hands me for the other three controls.

LEMMA (wrap budget).  m_h(x) = w - (x_1-x_h) satisfies m_h(0)=w, m_h>=0 on R, and
m_h(x+theta) - m_h(x) = theta_h - theta_1.  If J_open(theta)=empty then
theta = 1^a 0^{h-a}, so theta_h-theta_1 = -1 when 0<a<h and 0 when a in {0,h}.
Hence along an all-J_open-empty path,  #{a : 0 < alpha_a < h} <= w.

TEST 1: verify the budget inequality directly over many (h,w,alpha) -- INCLUDING
        parameter points where all-J_open-empty tableaux EXIST (firing control).
TEST 2: verify m_i >= 1 for every i in J_open(theta) on every legal strip (so the
        factor (1-t^{m_i}) genuinely vanishes at t=1 -- korff.py's assert, as a claim).
TEST 3: THE SCREEN.  At h=2 and alpha=(1^n), theta is a UNIT vector.  Any weight of
        index-set type dies at t=1 unless its index set is EMPTY on BOTH unit vectors.
        Evaluate all five variants of controls.py on (1,0) and (0,1).
TEST 4: general k -- is c_alpha(1) = A-B nonzero for k=2 too?  (the refutation needs it)
"""
import sys, itertools
sys.path.insert(0, '/home/clio/projects/proofs/code-1010-q410')
from korff import (paths, J_open, J_cyclic, I_cyclic, mvec, in_region,
                   step_vectors, weight_step, weight_path, t)
from controls import weight_step_variant
import sympy as sp

def thetas(P):
    h = len(P[0])
    return [tuple(P[a+1][i]-P[a][i] for i in range(h)) for a in range(len(P)-1)]

print("="*78)
print("TEST 1 -- the wrap-budget inequality  #{a : 0<alpha_a<h} <= w  on all-J_open-empty T")
print("="*78)
fire = 0; viol = 0; checked = 0
for h in [2,3,4]:
    for w in [2,3,4]:
        for n in [2,3,4]:
            for alpha in itertools.product(range(h+1), repeat=n):
                for P in paths(h, w, alpha):
                    if any(J_open(th) for th in thetas(P)): continue
                    checked += 1; fire += 1
                    nnc = sum(1 for a in alpha if 0 < a < h)
                    if nnc > w:
                        viol += 1
                        print(f"   VIOLATION h={h} w={w} alpha={alpha} path={P}")
print(f"   all-J_open-empty tableaux found (FIRING CONTROL, must be > 0): {fire}")
print(f"   budget-inequality violations among them: {viol}")
# negative control: the same count WITHOUT the J_open=empty filter must violate
nv = 0; tot = 0
for h in [2,3]:
    for w in [2,3]:
        for n in [4]:
            for alpha in itertools.product(range(h+1), repeat=n):
                for P in paths(h, w, alpha):
                    tot += 1
                    if sum(1 for a in alpha if 0 < a < h) > w: nv += 1
print(f"   NEGATIVE CONTROL: without the J_open=empty filter, {nv}/{tot} tableaux")
print(f"   violate the same inequality -- so the filter is doing the work, not the bound.")

print()
print("="*78)
print("TEST 2 -- m_i >= 1 for every i in J_open(theta) on every legal strip")
print("="*78)
bad = 0; seen = 0
for h in [2,3,4,5]:
    for w in [2,3,4,5]:
        reps = [x for x in itertools.product(range(w+1), repeat=h) if in_region(x,w)]
        for x in reps:
            for a in range(h+1):
                for th in step_vectors(h, a):
                    q = tuple(x[i]+th[i] for i in range(h))
                    if not in_region(q, w): continue
                    m = mvec(x, w)
                    for i in J_open(th):
                        seen += 1
                        if m[i] < 1: bad += 1; print(f"   BAD x={x} th={th} i={i} m={m}")
print(f"   (i, strip) pairs with i in J_open: {seen}   with m_i < 1: {bad}")
# firing control: the SAME claim for J_cyclic's wrap index i=h-1 -- does it also hold?
badc = 0; seenc = 0
for h in [2,3,4]:
    for w in [2,3,4]:
        for x in [x for x in itertools.product(range(w+1), repeat=h) if in_region(x,w)]:
            for a in range(h+1):
                for th in step_vectors(h, a):
                    q = tuple(x[i]+th[i] for i in range(h))
                    if not in_region(q, w): continue
                    m = mvec(x, w)
                    for i in J_cyclic(th):
                        seenc += 1
                        if m[i] < 1: badc += 1
print(f"   CONTROL on J_cyclic: {seenc} pairs, {badc} with m_i < 1 "
      f"(expected 0 too -- Korff's region does force it)")

print()
print("="*78)
print("TEST 3 -- THE SCREEN.  h=2, alpha=(1^n): theta is a unit vector.")
print("   An index-set weight survives t=1 only if its index set is EMPTY on BOTH")
print("   unit vectors.  Below: index set and Psi(1) on (1,0) and (0,1), at a state")
print("   x=(2,0) with w=4 where both steps are legal.")
print("="*78)
h = 2; w = 4; x = (2,0)
print(f"   state x={x}, w={w}, m(x)={mvec(x,w)};  legal: "
      f"{[th for th in [(1,0),(0,1)] if in_region(tuple(x[i]+th[i] for i in range(2)), w)]}")
rows = []
for variant in ["true_Psi","open_J","shift_m","use_mlam","I_with_mu"]:
    cells = []
    for th in [(1,0),(0,1)]:
        q = tuple(x[i]+th[i] for i in range(h))
        if not in_region(q, w): cells.append("illegal"); continue
        W = weight_step_variant(x, th, w, variant)
        cells.append(f"{W}  -> at t=1: {sp.expand(W).subs(t,1)}")
    survives = all("t=1: 0" not in c for c in cells)
    rows.append((variant, cells, survives))
    print(f"   {variant:10s} theta=(1,0): {cells[0]:34s} | theta=(0,1): {cells[1]}")
    print(f"   {'':10s}   -> survives the t=1 screen on BOTH unit steps? {survives}")
print(f"   VARIANTS SURVIVING: {[r[0] for r in rows if r[2]]}")

print()
print("="*78)
print("TEST 4 -- general k: is c_alpha(1) = A_alpha - B_alpha nonzero?")
print("   (the refutation needs c_alpha(1) != 0; k=1 is tabulated, k=2 is not)")
print("="*78)
def cminus(la, k, w):
    m = 2*k
    if all(la[2*i] == la[2*i+1] for i in range(k)): return 1
    if la[0]-la[m-1] == w and all(la[2*i+1] == la[2*i+2] for i in range(k-1)): return -1
    return 0
for (k, ell) in [(1,0),(1,1),(1,2),(2,0),(2,1),(3,0)]:
    h = 2*k; w = 2*ell+2; M = k+ell+1; n = 2*M
    alpha = (1,)*n
    A=B=0; nwd=0; psi1 = sp.Integer(0)
    for P in paths(h, w, alpha):
        if not any(J_open(th) for th in thetas(P)): nwd += 1
        c = cminus(P[-1], k, w)
        A += (c==1); B += (c==-1)
        if c: psi1 += c*sp.expand(weight_path(P, w, "Psi_open")).subs(t,1)
    print(f"   k={k} ell={ell} h={h} w={w} M={M} n={n}:  A={A} B={B}  c_alpha(1)=A-B={A-B}"
          f"   all-J_open-empty T: {nwd}   sum_T c^- Psi_open_T(1) = {psi1}"
          f"   {'REFUTED' if A-B != psi1 else 'consistent'}")
