from puzzles import *
from kalg import *
import itertools, sys
KTW   = ['zero_A','zero_B','one_A','one_B','rh_014','rh_230','rh_125']
SET_HEX = KTW + ['hex']
# section ss:protrude complete set, verbatim from the figure (ktiles.tex:1601-1625):
# 2 zero-triangles, 3 one-triangles, 4 rhombi, 3 green pieces = 12 tiles.
PROT = ['zero_A','zero_B','one_A','one_A2','one_B2',
        'rh_014','rh_230','rh_125a','rh_125b','hex','hex2','trap']
def parts_in(rows, cols):
    out=[]
    def rec(pref,prev,r):
        if r==0: out.append(tuple(pref)); return
        for v in range(min(prev,cols),-1,-1): rec(pref+[v],v,r-1)
    rec([],cols,rows); return [tuple(x for x in p if x>0) for p in out]
def pad(p,k): return list(p)+[0]*(k-len(p))

k=int(sys.argv[1]); m=int(sys.argv[2])     # ambient = k rows x m cols
n=k+m                                      # "natural" side length
P=parts_in(k,m)
print("### ambient %d rows x %d cols ; refined set at side n=%d  vs  \\hex at side n+1=%d" % (k,m,n,n+1))
dc={nu:d_coeffs(nu,max(2,k),max(2,k),maxdeg=2*max(1,sum(nu))) for nu in P}
rows=nz=b1=b2=0
for la,mu in itertools.product(P,P):
    for nu in P:
        if sum(nu)>sum(la)+sum(mu): continue
        exp=dc[nu].get((tuple(x for x in la if x>0),tuple(x for x in mu if x>0)),0)
        g1=count_puzzles(n,  pad(la,k),pad(mu,k),pad(nu,k),k,PROT)      # refined, side n
        g2=count_puzzles(n+1,pad(la,k),pad(mu,k),pad(nu,k),k,SET_HEX)   # \hex, side n+1
        rows+=1
        if exp!=0: nz+=1
        if exp!=g1: b1+=1; print("  PROT  MISMATCH la=%s mu=%s nu=%s  d=%d prot=%d"%(la,mu,nu,exp,g1))
        if exp!=g2: b2+=1; print("  HEX   MISMATCH la=%s mu=%s nu=%s  d=%d hex=%d"%(la,mu,nu,exp,g2))
print("  refined set (side %d) vs d : %d/%d agree"%(n,rows-b1,rows))
print("  \\hex       (side %d) vs d : %d/%d agree"%(n+1,rows-b2,rows))
print("  NONZERO d among those rows : %d/%d  (so %d rows are 0=0)"%(nz,rows,rows-nz))
