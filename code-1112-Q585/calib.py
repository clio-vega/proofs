import sys; sys.path.insert(0,'.'); sys.path.insert(0,'/home/clio/git/puzzles')
from puzzle import *; from commutor import parts, LR
from solve_puzzle import partition_to_binary
from itertools import product
def enc(p,n,m): return list(partition_to_binary(p,n,m))
def cases(n,m):
    for a in range(1,n*m+1):
      for lam in parts(a,n):
        if lam[0]>m: continue
        for b in range(1,n*m+1-a+1):
          for mu in parts(b,n):
            if mu[0]>m: continue
            L=LR(lam,mu,n)
            for nu in parts(a+b,n):
              if nu[0]>m: continue
              yield lam,mu,nu,len(L.get(tuple(list(nu)+[0]*(n-len(nu))),[]))
if __name__=='__main__':
    n,m=3,3; C=list(cases(n,m))
    for rS,rW,rE in product((0,1),repeat=3):
      for k,f in readings.items():
        ok=0
        for lam,mu,nu,t in C:
            S=enc(nu,n,m); W=enc(lam,n,m); E=enc(mu,n,m)
            if rS: S=S[::-1]
            if rW: W=W[::-1]
            if rE: E=E[::-1]
            ok+= len(puzzles(S,W,E,f))==t
        print(rS,rW,rE,k,ok,'/',len(C))
