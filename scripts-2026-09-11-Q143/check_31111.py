from matrix_engine import two_term_relations, matrix, is_rim_hook
from check_matrix import col2pair
from hookgraph import hooks
from graph import edges, vertices, partitions, maya, components_and_bipartite
from itertools import combinations

mu=(3,1,1,1,1)
print('hooks of',mu,':',{c:h for c,h in sorted(hooks(mu).items())})
cols,rels=two_term_relations(mu)
P={g:col2pair(mu,g) for g in cols}
inv={v:g for g,v in P.items()}
# cells (1,2),(1,3) <-> pairs (2,0),(2,1)  [b_1=2, holes -5,0,1]
tgt=frozenset(((2,0),(2,1)))
found=[(g1,g2,h) for (g1,g2,h) in rels if frozenset((P[g1],P[g2]))==tgt]
print('E2 (matrix, rim hooks only) two-term rows joining cells (1,2)-(1,3):',found)
print('   i.e. columns gamma =',[(g1,g2) for g1,g2,_ in found])
rows,colsM,M=matrix(mu)
for (g1,g2,h) in found:
    for lam in rows:
        nz=[(g,M[(lam,g)]) for g in colsM if (lam,g) in M]
        if len(nz)==2 and {nz[0][0],nz[1][0]}=={g1,g2}:
            print('   witnessing row lambda =',lam,' entries',nz)
print()
print('paper: "no second edge inside row 1 exists" for (3,1,1,1,1).')
print('  sufficient condition d1+d2 > h_11 :', 2+1,'>',7,'=',2+1>7)
print('  EXACT condition h_12+h_13 =',3,'is a hook of row 1?',3 in {7,2,1},'-> edge EXISTS')
print('  row-1 triangle {(1,1),(1,2),(1,3)} present in E1:',
      all(frozenset((x,y)) in {frozenset((a,b)) for a,b,_,_ in edges(mu)}
          for x,y in [((2,-5),(2,0)),((2,-5),(2,1)),((2,0),(2,1))]))
print()
# conjugation duality of the two controls
def conjp(mu):
    return tuple(sum(1 for x in mu if x>=j) for j in range(1,mu[0]+1))
bad=0
for n in range(1,10):
    for mu in partitions(n):
        V=vertices(mu); Vc=vertices(conjp(mu))
        a=components_and_bipartite(V,edges(mu,use_I=False))
        b=components_and_bipartite(Vc,edges(conjp(mu),use_II=False))
        if (len(a[0]),bool(a[3]))!=(len(b[0]),bool(b[3])): bad+=1; print('dual fail',mu)
print('CTRL-1(mu) vs CTRL-2(mu\'): (#components, non-bipartite) agree for all mu, n<=9; failures =',bad)
print('  -> the two controls have equal TOTALS by conjugation, not by being the same test.')
print('  per-mu they differ: e.g. mu=(3,1): CTRL-1 #comps =',
      len(components_and_bipartite(vertices((3,1)),edges((3,1),use_I=False))[0]),
      ', CTRL-2 #comps =',len(components_and_bipartite(vertices((3,1)),edges((3,1),use_II=False))[0]))
