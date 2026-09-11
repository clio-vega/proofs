from graph import *
from fractions import Fraction
import sympy as sp
t=sp.symbols('t')

def cycle_space_labels(V,E):
    """Return True if SOME cycle has label product != 1.  Spanning forest + fundamental cycles
    generate the cycle space; product is multiplicative, so 'all fundamental cycles trivial'
    <=> 'all cycles trivial'.  Label of edge (x,y,h): w(y) = -t^h w(x)."""
    adj={v:[] for v in V}
    for i,(x,y,h,ty) in enumerate(E):
        adj[x].append((y,i,+1)); adj[y].append((x,i,-1))
    par={}; seen=set(); tree=set()
    # potential: phi(v) = (parity, exponent) so that w(v) = (-1)^par t^exp w(root)
    pot={}
    for s in V:
        if s in seen: continue
        seen.add(s); pot[s]=(0,0); stack=[s]
        while stack:
            x=stack.pop()
            for (y,i,sg) in adj[x]:
                if y not in seen:
                    seen.add(y); tree.add(i)
                    p,e=pot[x]
                    pot[y]=((p+1)%2, e+sg*E[i][2])
                    stack.append(y)
    for i,(x,y,h,ty) in enumerate(E):
        if i in tree: continue
        px,ex=pot[x]; py,ey=pot[y]
        # cycle label = (-1)^{1+px+py} t^{h + ex - ey}  ... check: w(y)=-t^h w(x) vs tree path
        par=(1+px+py)%2; exp=h+ex-ey
        if par==1 or exp!=0:
            return True
    return False

rows={}
EXC=[]
for n in range(1,10):
    ps=list(partitions(n)); nc=0; nb=0; ncrit=0; bad=[]
    for mu in ps:
        V=vertices(mu); E=edges(mu)
        assert len(V)==n, (mu,len(V))
        comps,comp,colour,odd=components_and_bipartite(V,E)
        conn = (len(comps)==1)
        nonbip = (len(odd)>0)
        crit = cycle_space_labels(V,E)
        if conn: nc+=1
        if nonbip: nb+=1
        if crit: ncrit+=1
        else: bad.append(mu)
        # consistency: nonbipartite => criterion holds
        assert (not nonbip) or crit, mu
    rows[n]=(len(ps),nc,nb,ncrit,bad)
    print(f"n={n:2d}  #mu={len(ps):3d}  connected={nc:3d}  non-bipartite={nb:3d}  cycle-crit={ncrit:3d}  FAIL={bad}")
print()
print("totals:", sum(r[0] for r in rows.values()), sum(r[1] for r in rows.values()),
      sum(r[2] for r in rows.values()), sum(r[3] for r in rows.values()))
