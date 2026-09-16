"""Is R_e(t) a BOSON?  KMS/Uglov B_{-1} satisfies [B_1,B_{-1}] = scalar (a q-integer).
Test the analogous statement for R_e with the Hall pairing:  R_e^* = "remove a
connected e-ribbon, weight t^ht".  Is [R_e^*, R_e] a scalar operator?

Also: verify (1+t) | [R_e,R_f], which C1 predicts (R_e(-1)=M_{p_e} all commute)."""
import sys, sympy as sp
sys.path.insert(0, '/home/clio/projects/probes/2026-09-06-Q84')
import engine
from bead import trim, t

def Rstar(state, e):
    """adjoint of R_e w.r.t. <s_lam,s_mu>=delta: remove a connected e-ribbon."""
    out = {}
    for mu, w in state.items():
        for lam in engine.parts_of(sum(mu) - e) if sum(mu) >= e else ():
            if not engine.contains(mu, lam): continue
            cells = engine.boxes(mu, lam)
            if len(cells) != e: continue
            m, ht, sq = engine.components_and_ht(cells)
            if sq or m != 1: continue
            out[lam] = out.get(lam, 0) + w * t**ht
    return engine.clean(out)

print("=== [R_e^*, R_e] s_lam  (scalar <=> boson) ===")
for e in (2, 3):
    for n in range(0, 6):
        for lam in (engine.parts_of(n) if n else ((),)):
            a = Rstar(engine.op_R({trim(lam): 1}, e), e)
            b = engine.op_R(Rstar({trim(lam): 1}, e), e)
            c = engine.sub(a, b)
            diag = sp.factor(c.get(trim(lam), 0))
            offd = {k: sp.factor(v) for k, v in c.items() if k != trim(lam)}
            if offd or n <= 2:
                print(f"  e={e} lam={lam}: diag={diag}  off-diagonal={offd if offd else '{}'}")
        if n >= 3: break
    print()

print("=== (1+t) | [R_e,R_f] ? ===")
bad = 0; tot = 0
for (e, f) in [(1,2),(1,3),(2,3),(2,4),(3,4),(3,5)]:
    for n in range(0, 6):
        for lam in (engine.parts_of(n) if n else ((),)):
            A = engine.op_R(engine.op_R({trim(lam): 1}, f), e)
            B = engine.op_R(engine.op_R({trim(lam): 1}, e), f)
            for k, v in engine.sub(A, B).items():
                tot += 1
                if sp.simplify(sp.expand(v).subs(t, -1)) != 0: bad += 1
print(f"  commutator coefficients tested: {tot};  nonzero at t=-1: {bad}")
