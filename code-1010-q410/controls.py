"""DISCRIMINATION CONTROLS.  Each 'variant' is a deliberately wrong reading of Korff
(5.12)/(5.13).  A gate is only informative if the variants FAIL it."""
import sys, itertools
sys.path.insert(0,'.')
from korff import *
import sympy as sp
t = sp.Symbol('t')

def weight_step_variant(x, theta, w, variant):
    h = len(x); m = mvec(x, w)
    mlam = mvec(tuple(x[i]+theta[i] for i in range(h)), w)
    if variant == "true_Psi":      idx, mm = J_cyclic(theta), m
    elif variant == "open_J":      idx, mm = J_open(theta), m          # no wrap
    elif variant == "shift_m":     idx, mm = J_cyclic(theta), tuple(m[(i+1)%h] for i in range(h))
    elif variant == "use_mlam":    idx, mm = J_cyclic(theta), mlam     # outer instead of inner
    elif variant == "I_with_mu":   idx, mm = I_cyclic(theta), m        # wrong index set
    else: raise ValueError(variant)
    out = sp.Integer(1)
    for i in idx:
        if mm[i] < 1: return sp.Integer(0)      # variants may hit m=0; (1-t^0)=0
        out *= (1-t**mm[i])
    return sp.expand(out)

def b(m):
    out = sp.Integer(1)
    for mi in m:
        for j in range(1, mi+1): out *= (1-t**j)
    return sp.expand(out)

def gate_lemma54_variant(h, w, variant):
    """Korff eq (5.16):  b_lam * Psi_{lam/d/mu} = b_mu * Phi_{lam/d/mu}.
    Phi is held at Korff's true (5.12); only Psi is varied."""
    checked = fail = 0
    reps = [x for x in itertools.product(range(w+1), repeat=h) if x[h-1]==0 and in_region(x,w)]
    for x in reps:
        for a in range(h+1):
            for th in step_vectors(h, a):
                q = tuple(x[i]+th[i] for i in range(h))
                if not in_region(q, w): continue
                Ps = weight_step_variant(x, th, w, variant)
                Ph = weight_step(x, th, w, "Phi")
                checked += 1
                if sp.expand(b(mvec(q,w))*Ps - b(mvec(x,w))*Ph) != 0: fail += 1
    return checked, fail

def gate_mac_order_variant(h, nvars, wbig, variant):
    """Order at which the variant's weight departs from Macdonald's ordinary psi.
    The TRUE cylindric weight must depart only through the wrap factor (1-t^{m_h}),
    m_h = w-(x_1-x_h), so the departure order is >= w - span: LARGE.
    A misplaced wrap departs at LOW order."""
    worst = None
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, wbig, alpha):
            A = sp.Integer(1); B = weight_path(P, wbig, "Psi_open")
            for a in range(len(P)-1):
                th = tuple(P[a+1][i]-P[a][i] for i in range(h))
                A *= weight_step_variant(P[a], th, wbig, variant)
            d = sp.expand(A-B)
            if d != 0:
                lo = min(mm[0] for mm in sp.Poly(d,t).monoms())
                worst = lo if worst is None or lo < worst else worst
    return worst

if __name__ == "__main__":
    print("GATE 2 (Korff eq (5.16), b_lam*Psi = b_mu*Phi) WITH CONTROLS")
    for variant in ["true_Psi","open_J","shift_m","use_mlam","I_with_mu"]:
        tot = bad = 0
        for (h,w) in [(2,2),(3,2),(2,3),(3,3),(4,2)]:
            c,f = gate_lemma54_variant(h,w,variant); tot += c; bad += f
        flag = "PASS" if bad==0 else "fail"
        print(f"   {variant:12s}: {bad:4d}/{tot} strips violate (5.16)   [{flag}]")
    print()
    print("GATE 3 (departure order from Macdonald's ordinary psi; large = wrap only) WITH CONTROLS")
    for variant in ["true_Psi","shift_m","use_mlam","I_with_mu"]:
        o1 = gate_mac_order_variant(2,3,8,variant)
        o2 = gate_mac_order_variant(3,3,9,variant)
        print(f"   {variant:12s}: lowest departure degree  h=2,w=8: {o1}   h=3,w=9: {o2}")
