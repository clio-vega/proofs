from graph import partitions, vertices, edges, components_and_bipartite
from hookgraph import hooks, conj, hook_edges
from itertools import combinations

def line_tri(hs):
    S=set(hs); return any(hs[q]+hs[s] not in S for q,s in combinations(range(1,len(hs)),2))
def graph_triangle_free(mu):
    h=hooks(mu); mp=conj(mu)
    L=[[h[(i,j)] for j in range(1,mu[i-1]+1)] for i in range(1,len(mu)+1)]
    L+=[[h[(i,k)] for i in range(1,mp[k-1]+1)] for k in range(1,mu[0]+1)]
    return not any(line_tri(x) for x in L)
def real_triangle_free(mu):
    """brute force over all 3-subsets of cells -- does NOT use the line lemma"""
    E={frozenset((x,y)) for (x,y,ty) in hook_edges(mu)}
    V=[(i,k) for i in range(1,len(mu)+1) for k in range(1,mu[i-1]+1)]
    return not any(frozenset((a,b)) in E and frozenset((a,c)) in E and frozenset((b,c)) in E
                   for a,b,c in combinations(V,3))
EXC={(2,2),(3,2),(2,2,1)}
def predict(mu):   # Theorem C
    return sum(mu)<=3 or mu in EXC

bad=0; nb=0; tot=0
for n in range(1,13):
    for mu in partitions(n):
        tot+=1
        t1=graph_triangle_free(mu); t2=real_triangle_free(mu); p=predict(mu)
        if not (t1==t2==p): bad+=1; print('THM-C FAIL',mu,t1,t2,p)
        if n<=10:
            V=vertices(mu); bip = not components_and_bipartite(V,edges(mu))[3]
            if bip!=t2: nb+=1; print('bipartite != triangle-free',mu)
print(f'Theorem C: line-lemma == brute-force-triangle == prediction, on {tot} partitions n<=12; failures={bad}')
print(f'bipartite <=> triangle-free on all mu, n<=10; failures={nb}')

# the n=2*ell+1 characterisation of row-1 triangle-freeness
f=0;c=0
for n in range(1,13):
    for mu in partitions(n):
        h=hooks(mu); ell=len(mu)
        r1=[h[(1,j)] for j in range(1,mu[0]+1)]
        pred = (mu[0]<=2) or (mu[0]==3 and n==2*ell+1)
        c+=1
        if (not line_tri(r1))!=pred: f+=1; print('row1 char fail',mu)
print(f'row 1 triangle-free <=> mu_1<=2 or (mu_1=3 and n=2*ell+1): {c-f}/{c} n<=12')

# signed exponent sums on the three exceptions (boundary population, explicit)
import sympy as sp
t=sp.symbols('t')
for mu in [(2,2),(3,2),(2,2,1),(3,),(2,1),(1,1,1)]:
    V=vertices(mu); E=edges(mu)
    comps,comp,colour,odd=components_and_bipartite(V,E)
    # fundamental cycles
    adj={v:[] for v in V}
    for i,(x,y,hh,ty) in enumerate(E): adj[x].append((y,i,+1)); adj[y].append((x,i,-1))
    pot={}; seen=set(); tree=set()
    for s in V:
        if s in seen: continue
        seen.add(s); pot[s]=(0,0); st=[s]
        while st:
            x=st.pop()
            for (y,i,sg) in adj[x]:
                if y not in seen:
                    seen.add(y); tree.add(i); p,e=pot[x]; pot[y]=((p+1)%2,e+sg*E[i][2]); st.append(y)
    cyc=[]
    for i,(x,y,hh,ty) in enumerate(E):
        if i in tree: continue
        px,ex=pot[x]; py,ey=pot[y]
        cyc.append(((1+px+py)%2, hh+ex-ey))
    print(f'{str(mu):12s} |V|={len(V)} |E|={len(E)} bipartite={not odd} '
          f'fundamental cycles (parity, signed exp sum) = {cyc}')
