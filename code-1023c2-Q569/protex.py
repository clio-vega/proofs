from puzzles import *
from kalg import *
KTW=['zero_A','zero_B','one_A','one_B','rh_014','rh_230','rh_125']
PROT=['zero_A','zero_B','one_A','one_A2','one_B2','rh_014','rh_230','rh_125a','rh_125b','hex','hex2','trap']
la,mu,nu,k=[2,1,0],[4,2,0],[4,3,1],3
dc=d_coeffs(tuple(nu),3,3,maxdeg=2*sum(nu))
d=dc.get(((2,1),(4,2)),0)
print("paper's ss:protrude example  la=(2,1) mu=(4,2) nu=(4,3,1), k=3")
print("  independent algebra: d^nu_{la,mu} =",d)
for n in (7,8):
    print("  \\hex-puzzles at side n=%d : %d   (thm:hex-tiling needs partitions in %d rows x %d cols; la_1,mu_1,nu_1<=%d)"
          % (n, count_puzzles(n,la,mu,nu,k,KTW+['hex']), k, n-k-1, n-k-1))
print("  refined (ss:protrude) set at side 7 : %d" % count_puzzles(7,la,mu,nu,k,PROT))
