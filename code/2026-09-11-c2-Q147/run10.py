import sympy as sp
from run6 import sector_eqs
for (e,f,mx) in [(2,3,11),(2,4,11),(3,4,12),(3,5,13)]:
    one, two = sector_eqs(e,f,mx, shared=False)
    syms = sorted({s for q in one|two for s in q.free_symbols}, key=str)
    print("e,f=(%d,%d) INDEPENDENT weights, syms=%s" % (e,f,syms))
    print("   1-bead locus:", sp.solve(list(one), syms, dict=True))
    print("   2-bead locus:", sp.solve(list(two), syms, dict=True))
