"""Q473 main verification: the residue functional R_kappa and its triangularity.

R_kappa = Res_{x_d = 1/(q y_{kappa_d})} ... Res_{x_1 = 1/(q y_{kappa_1})},  d = |kappa|.

PREDICTED (proved structurally):  R_kappa[G_{kappa'}] = 0  unless  kappa' >= kappa,
where kappa' >= kappa  iff  |kappa'| >= |kappa| and kappa'_i >= kappa_i for i <= |kappa|;
and R_kappa[G_kappa] != 0.
"""
import engine as E, sympy as sp, random, sys
from sympy import Rational as R
q, a = E.q, E.a

def numeric_env(L, maxy, seed=0):
    rr = random.Random(seed)
    s = {q: R(5, 97), a: R(37, 101)}
    for j in range(1, maxy+1):
        s[E.ysym(j)] = R(rr.randint(3, 29), 103+j)   # y's distinct, generic
    return s

def R_functional(expr, kappa, env):
    """iterated residue in x_1..x_d at x_i = 1/(q y_{kappa_i})"""
    f = expr
    for i, m in enumerate(kappa, start=1):
        xi = E.xsym(i)
        xi0 = 1/(env[q]*env[E.ysym(m)])
        f = sp.cancel(sp.together((xi - xi0)*f))
        f = sp.cancel(f.subs(xi, xi0))
    return sp.cancel(f)

def run(L, maxel, maxcard, seed=0):
    env = numeric_env(L, maxel, seed)
    kaps = [k for k in E.subsets_upto(maxel, maxcard) if len(k) <= L]
    # free variables x_{d+1}..x_L get generic rational values AFTER the residues;
    # assign them up front (they are not residue variables for any kappa of max card)
    rr = random.Random(seed+500)
    free = {E.xsym(i): R(rr.randint(7, 61), 107+i) for i in range(1, L+1)}
    Gs = {}
    for kp in kaps:
        Gs[kp] = sp.cancel(E.G(kp, L).subs(env))
    rows = []
    bad = []
    for kap in kaps:
        if len(kap) == 0:
            continue
        d = len(kap)
        for kp in kaps:
            sub = {E.xsym(i): free[E.xsym(i)] for i in range(d+1, L+1)}
            val = R_functional(Gs[kp].subs(sub), kap, env)
            iszero = (sp.simplify(val) == 0)
            pred_nonzero = E.dominates(kp, kap)
            agree = (iszero != pred_nonzero)
            rows.append((kap, kp, 0 if iszero else 1, 1 if pred_nonzero else 0, agree))
            if not agree:
                bad.append((kap, kp, iszero, pred_nonzero))
        sys.stdout.flush()
    ok = sum(1 for r in rows if r[4])
    print("L=%d maxel=%d maxcard=%d : %d/%d cells match prediction" % (L, maxel, maxcard, ok, len(rows)))
    nz = sum(1 for r in rows if r[2] == 1); z = len(rows)-nz
    print("   zero cells %d, nonzero cells %d  (so the prediction is NOT vacuous either way)" % (z, nz))
    diag = [(k, sp.simplify(R_functional(Gs[k].subs({E.xsym(i): free[E.xsym(i)] for i in range(len(k)+1, L+1)}), k, env)) != 0) for k in kaps if len(k) > 0]
    print("   diagonal nonzero: %d/%d" % (sum(1 for _, t in diag if t), len(diag)))
    for b in bad[:10]:
        print("   MISMATCH", b)
    return rows, bad

if __name__ == '__main__':
    L = int(sys.argv[1]); maxel = int(sys.argv[2]); maxcard = int(sys.argv[3])
    run(L, maxel, maxcard)
