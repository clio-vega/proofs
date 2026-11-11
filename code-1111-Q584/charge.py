"""Q584 C2: conserved 'charge' omega(d,label) in Q^m with omega(d+3,l) = -omega(d,l).
Every allowed tile (ccw edge list) must have sum omega = 0  (mod a constant c per triangle orientation
is allowed too: the sl3 invariant has weight (1,1,1); we model that by allowing
sum_up = c, sum_dn = -c for an unknown vector c ... handled by an extra unknown 'c').
Solve the linear system over Q for scalar omega (vector omega = independent copies), get the nullspace,
then evaluate the rim sums of hex, hexR, trap(+wall), K, tri on every basis element."""
import sys, itertools
sys.path.insert(0,__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','code-1023c2-Q569'))
from puzzles import TILES
from sympy import Matrix, Rational
X='X'; LABS=(0,1,X)
# unknowns: omega(d,l) for d in 0,1,2 (d>=3 is minus), plus c
idx={(d,l):i for i,(d,l) in enumerate(itertools.product(range(3),LABS))}; NC=len(idx)
def row(edges, cs=0):
    r=[0]*(NC+1)
    for d,l in edges:
        s=1 if d<3 else -1; r[idx[(d%3,l)]]+=s
    r[NC]-=cs        # sum omega - cs*c = 0
    return r
def cyc(t): return {t,t[1:]+t[:1],t[2:]+t[:2]}
UP=set().union(cyc((0,0,0)),cyc((1,1,1)),cyc((1,X,0)),{(X,X,X)})
DN=set().union(cyc((0,0,0)),cyc((1,1,1)),cyc((0,1,X)),{(X,X,X)})
rows=[]
for t in UP: rows.append(row(zip((0,2,4),t),1))
for t in DN: rows.append(row(zip((1,3,5),t),-1))
M=Matrix(rows); N=M.nullspace()
print('triangles:',len(rows),' unknowns:',NC+1,' nullspace dim:',len(N))
def rim(nm): return row(TILES[nm],0)
names=['zero_A','one_A','rh_014','rh_230','rh_125','K','tri','hex','hexR']
for nm in names:
    vals=[sum(a*b for a,b in zip(rim(nm),v)) for v in N]
    print(nm.ljust(7),'area',None,' rim-sum on nullspace basis:',vals)
# describe nullspace
for v in N: print('basis',{k:v[i] for k,i in idx.items() if v[i]!=0},'c=',v[NC])
