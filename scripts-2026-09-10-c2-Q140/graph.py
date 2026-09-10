"""The relation graph G on the n columns (= cells of mu).
Each edge:  w(v) = -t^h w(u).  Assign a "potential" phi(v) in Z x Z/2 (t-exponent, sign).
w = 0 is forced  <=>  every component carries a cycle with nontrivial monodromy.
"""
from rim import *
from reflect import derived_relations, cols_abacus

def analyse(mu):
    n = sum(mu); N = n + len(mu) + 2
    B, cabac, rels = derived_relations(mu, N)
    V = list(cabac.keys())
    idx = {v:i for i,v in enumerate(V)}
    adj = {v: [] for v in V}
    for kind, u, v, h in rels:
        adj[u].append((v, h, +1)); adj[v].append((u, -h, -1))   # w(v)=-t^h w(u)
    # BFS assigning potential (exponent, parity); detect nontrivial monodromy
    seen = {}; comps = []
    for s in V:
        if s in seen: continue
        comp = [s]; seen[s] = (0,0); dead = False; stack=[s]
        while stack:
            u = stack.pop()
            for (v,h,sg) in adj[u]:
                pot = (seen[u][0] + h, (seen[u][1] + 1) % 2)
                if v not in seen:
                    seen[v] = pot; comp.append(v); stack.append(v)
                elif seen[v] != pot:
                    dead = True
        comps.append((comp, dead))
    return V, comps

print(" n | mu                  | #comps | all comps dead? | rank(2-term)==n?")
for n in range(4, 10):
    for mu in partitions(n):
        V, comps = analyse(mu)
        alld = all(d for c,d in comps)
        rows, cols, M, b = local_system(mu)
        rr = [i for i,l in enumerate(rows)
              if l != tuple(mu) and sum(1 for j in range(M.cols) if M[i,j]!=0)==2]
        full = (M[rr,:].rank() == n) if rr else False
        mark = "" if alld == full else "   *** MISMATCH ***"
        print(f"{n:2d} | {str(mu):19s} | {len(comps):2d} | {str(alld):5s} | {str(full):5s}{mark}")
