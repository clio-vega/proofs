"""
Q460 control 1 (SATISFIABILITY, sec 6.1): reproduce the paper's three printed
N=2 expansions (Appendix chapter:examples) term by term from my own classifier.
Control 2 (NEGATIVE, sec 6.5): a planted swap of cases B<->C must FAIL.

Built from: the .tex source of 2412.18085 only.  The classifier in ikdelta.py
was written from eq:IK-delta1/2/3; the targets below are typed from the
APPENDIX expansions, which are a different part of the document.  So the two
sides are not built from each other.
"""
import sympy as sp, sys
sys.path.insert(0, '/home/clio/projects/code/q460')
import ikdelta as ik
from ikdelta import q, Delta, row_to_endpoints

x1, x2 = sp.symbols('x1 x2')
z1, z2 = ik.Z(1), ik.Z(2)

# ---- printed targets, typed from the appendix expansions ------------------
# S=(2,2), eq:=example--(2,2)
T22 = [
    (((2,0),(0,2)), (x2-q**2*x1)*(x2-q**3*x1)/((x2-x1)*(x2-q*x1))),
    (((1,1),(1,1)), (x1*x2+q**2*x1*x2-q**3*x2*z1-q**3*x1*z1)
                   *(x1*x2+q**2*x1*x2-q**2*x2*z2-q**2*x1*z2)
                   /(q*(1-q)**2*(x2-q*x1)*(x1-q*x2)*z1*z2)),
    (((0,2),(2,0)), (x1-q**2*x2)*(x1-q**3*x2)/((x1-x2)*(x1-q*x2))),
]
# S=(1,2,1), eq:=example--(1,2,1)
T121 = [
    (((1,0,1),(0,2,0)), (q**2*x2-x1)*(x2-q**2*x1)/(x2-x1)**2),
    (((1,1,0),(0,1,1)), (x2-q**2*x1)*(x1*x2+q**2*x1*x2-q**2*x2*z2-q**3*x1*z2)
                        /(q*(1-q)*(x2-x1)**2*z2)),
    (((0,1,1),(1,1,0)), (x1-q**2*x2)*(x1*x2+q**2*x1*x2-q**2*x1*z2-q**3*x2*z2)
                        /(q*(1-q)*(x1-x2)**2*z2)),
    (((0,2,0),(1,0,1)), (q**2*x1-x2)*(x1-q**2*x2)/(x1-x2)**2),
]
# S=(1,1,1,1), eq:=sym--(1,1,1,1)
T1111 = [
    (((1,1,0,0),(0,0,1,1)), (x2-q**2*x1)*(x2-q**3*x1)/((x2-x1)*(x2-q*x1))),
    (((1,0,1,0),(0,1,0,1)), (x2-q**2*x1)**2*(x1-q*x2)/((x2-x1)**2*(x2-q*x1))),
    (((1,0,0,1),(0,1,1,0)), (x2-q**2*x1)*(q**2*x2-x1)/(x2-x1)**2),
    (((0,1,1,0),(1,0,0,1)), (x1-q**2*x2)*(q**2*x1-x2)/(x1-x2)**2),
    (((0,1,0,1),(1,0,1,0)), (x1-q**2*x2)**2*(x2-q*x1)/((x1-x2)**2*(x1-q*x2))),
    (((0,0,1,1),(1,1,0,0)), (x1-q**2*x2)*(x1-q**3*x2)/((x1-x2)*(x1-q*x2))),
]

def run(label, targets, sym_used):
    ok = fail = 0
    for (U, V), target in targets:
        u, v = row_to_endpoints(list(U)), row_to_endpoints(list(V))
        mine = Delta(u, v, x1, x2, used_sym=sym_used)
        d = sp.simplify(sp.together(mine - target))
        if d == 0:
            ok += 1
        else:
            fail += 1
            print("   FAIL %s  U=%s V=%s  u=%s v=%s  residual=%s"
                  % (label, U, V, u, v, sp.srepr(d)[:120]))
    print("   %-12s %d/%d match" % (label, ok, ok+fail))
    return ok, fail

print("CONTROL 1 -- satisfiability: printed appendix expansions")
sym = set()
tot_ok = tot = 0
for lab, T in (("S=(2,2)", T22), ("S=(1,2,1)", T121), ("S=(1,1,1,1)", T1111)):
    o, f = run(lab, T, sym)
    tot_ok += o; tot += o+f
print("   TOTAL %d/%d" % (tot_ok, tot))
print("   cases exercised: %s" % sorted(sym))

# ---- CONTROL 2: planted swap B <-> C must break the match ----------------
print("\nCONTROL 2 -- negative control: plant a B<->C swap, it MUST fail")
ik.FN['B'], ik.FN['C'] = ik.FN['C'], ik.FN['B']
sym2 = set(); bad_ok = bad_tot = 0
for lab, T in (("S=(2,2)", T22), ("S=(1,2,1)", T121), ("S=(1,1,1,1)", T1111)):
    o, f = run(lab + "[planted]", T, sym2)
    bad_ok += o; bad_tot += o+f
print("   TOTAL under plant: %d/%d" % (bad_ok, bad_tot))
print("   NEGATIVE CONTROL %s" % ("FIRED (plant detected)" if bad_ok < tot_ok
                                  else "DID NOT FIRE -- controls are vacuous"))
