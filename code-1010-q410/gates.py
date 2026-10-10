"""Gates on the Korff-weight implementation.  Run BEFORE trusting any mismatch."""
import sys, itertools
sys.path.insert(0,'.')
from korff import *
import sympy as sp
t = sp.Symbol('t')

def gate_symmetry(h, w, nvars, verbose=True):
    """Korff's P_{lam/d/mu} is a symmetric function.  Fixed endpoint = fixed Korff shape.
    If my index set / wrap is wrong this fails."""
    bad = []; good = 0
    ends = {}
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, w, alpha):
            ends.setdefault(P[-1], []).append((alpha, P))
    xs = sp.symbols(f'x1:{nvars+1}', positive=True)
    for end, lst in sorted(ends.items()):
        F = sp.Integer(0)
        for alpha, P in lst:
            F += weight_path(P, w) * sp.prod([xs[i]**alpha[i] for i in range(nvars)])
        F = sp.expand(F)
        # test transposition x1<->x2 ... x_{nvars-1}<->x_{nvars}
        ok = True
        for i in range(nvars-1):
            sub = {xs[i]: xs[i+1], xs[i+1]: xs[i]}
            if sp.expand(F.subs(sub, simultaneous=True) - F) != 0:
                ok = False; break
        if ok: good += 1
        else: bad.append(end)
    if verbose:
        print(f"  gate_symmetry h={h} w={w} nvars={nvars}: {good} shapes symmetric, "
              f"{len(bad)} NOT symmetric {bad[:4]}")
    return len(bad) == 0, good, bad

def gate_symmetry_variant(h, w, nvars, kind):
    """Same test for a deliberately wrong weight -> must FAIL (discrimination control)."""
    ends = {}
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, w, alpha):
            ends.setdefault(P[-1], []).append((alpha, P))
    xs = sp.symbols(f'x1:{nvars+1}', positive=True)
    bad = 0
    for end, lst in sorted(ends.items()):
        F = sp.expand(sum(weight_path(P, w, kind) * sp.prod([xs[i]**a[i] for i in range(nvars)])
                          for a, P in lst))
        for i in range(nvars-1):
            sub = {xs[i]: xs[i+1], xs[i+1]: xs[i]}
            if sp.expand(F.subs(sub, simultaneous=True) - F) != 0:
                bad += 1; break
    print(f"  CONTROL kind={kind} h={h} w={w} nvars={nvars}: {bad} shapes NOT symmetric "
          f"(must be > 0 for a wrong weight)")
    return bad

def gate_lemma54(h, w):
    """Korff Lemma 5.4: b_lam/b_mu = Phi_{lam/d/mu} / Psi_{lam/d/mu},
    b_lam = prod_i (t)_{m_i(lam)},  (t)_m = (1-t)(1-t^2)...(1-t^m)."""
    def b(m):
        out = sp.Integer(1)
        for mi in m:
            for j in range(1, mi+1):
                out *= (1-t**j)
        return sp.expand(out)
    checked = fail = 0
    # all states x with x_h = 0 (fundamental representatives) and all single strips
    reps = [x for x in itertools.product(range(w+1), repeat=h)
            if x[h-1] == 0 and in_region(x, w)]
    for x in reps:
        for a in range(h+1):
            for th in step_vectors(h, a):
                q = tuple(x[i]+th[i] for i in range(h))
                if not in_region(q, w): continue
                Ps = weight_step(x, th, w, "Psi"); Ph = weight_step(x, th, w, "Phi")
                lhs = sp.expand(b(mvec(q,w))*Ps); rhs = sp.expand(b(mvec(x,w))*Ph)
                checked += 1
                if sp.simplify(lhs-rhs) != 0: fail += 1
    print(f"  gate_lemma54 h={h} w={w}: {checked} strips, {fail} failures of b_lam*Psi = b_mu*Phi")
    return fail == 0, checked

def gate_macdonald_limit(h, nvars, wbig):
    """Non-cylindric limit: for w >= (max column span)+1 the only difference between
    Korff's Psi and Macdonald's psi is the single wrap factor (1 - t^{m_h}),
    m_h = w - (x_1-x_h).  Check Psi = psi_open * (that factor when present)
    and that Psi == psi_open mod t^{w-span}, i.e. they agree to high order."""
    worst = None; cnt = 0
    for alpha in itertools.product(range(h+1), repeat=nvars):
        for P in paths(h, wbig, alpha):
            A = weight_path(P, wbig, "Psi"); B = weight_path(P, wbig, "Psi_open")
            cnt += 1
            d = sp.Poly(sp.expand(A-B), t)
            if d.total_degree() >= 0 and sp.expand(A-B) != 0:
                lo = min(m[0] for m in d.monoms() if d.coeff_monomial(m) != 0)
                worst = lo if worst is None or lo < worst else worst
    print(f"  gate_macdonald_limit h={h} nvars={nvars} w={wbig}: {cnt} tableaux; "
          f"Psi - psi_Macdonald has lowest surviving degree {worst} "
          f"(must be large: wrap factor only, order >= w-span)")
    return worst

if __name__ == "__main__":
    print("GATE 1  symmetry of Korff's P at fixed cylindric shape")
    for (h,w,nv) in [(2,2,3),(2,4,3),(3,2,3),(4,2,3),(3,3,3),(2,2,4),(4,2,4)]:
        gate_symmetry(h,w,nv)
    print("GATE 1c discrimination controls (wrong weights must break symmetry)")
    for kind in ["Psi_open","Phi"]:
        gate_symmetry_variant(3,2,3,kind); gate_symmetry_variant(4,2,3,kind)
    print("GATE 2  Korff Lemma 5.4")
    for (h,w) in [(2,2),(3,2),(2,3),(3,3),(4,2)]:
        gate_lemma54(h,w)
    print("GATE 3  non-cylindric limit vs Macdonald psi")
    gate_macdonald_limit(2,3,8); gate_macdonald_limit(3,3,9)
