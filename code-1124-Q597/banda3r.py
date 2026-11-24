"""VARIANT of banda3b.py with rhombi INCLUDED. Generalised reach of Theorem band for (A3): an A-shaped wake tile u read at k >= p+2 is COVERED if
 - u is a Lambda_k-triangle and proj(u), proj(u - alpha) lie in the same height group {t12,t23,t34} or {t45,t50,t01}
   at a band vertex x (Theorem band: k0=k2=k3=k5=0 there, so heights agree inside a group), or
 - u is a square over E_1(x) or E_4(x) and u - alpha is a square over the other (k1 = k4, same end heights).
usage: python3 banda3b.py SIZE"""
import sys
SIZE = sys.argv[1]
src = open('/home/clio/projects/proofs/code-1124-Q597/banda3.py').read().replace("if u not in W or u[0] == 'rh': continue", "if u not in W: continue")
src = src.replace("""            cov = False
            if zisL(u, k):
                pu = zpk(u, k)
                for x in bvs:
                    if pu == tri_at(x, want[0], 1): cov = True""", """            cov = False
            if zisL(u, k):
                pu = zpk(u, k); pa = zpk(ua, k)
                for x in bvs:
                    G1 = {tri_at(x, 1, 1), tri_at(x, 2, 1), tri_at(x, 3, 1)}; G2 = {tri_at(x, 4, 1), tri_at(x, 5, 1), tri_at(x, 0, 1)}
                    if (pu in G1 and pa in G1) or (pu in G2 and pa in G2): cov = 'same group'
            elif u[0] == 'sq':
                pu = zpk(u, k); pa = zpk(ua, k)
                for x in bvs:
                    E1 = frozenset((x, (x[0], x[1] + 1))); E4 = frozenset((x, (x[0], x[1] - 1)))
                    if {pu, pa} == {E1, E4}: cov = 'axis stacks'""")
assert "same group" in src and "if u not in W: continue" in src
exec(src)
