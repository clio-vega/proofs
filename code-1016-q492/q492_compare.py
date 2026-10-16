"""Do the two criteria agree on the 225 triples Q484 measured?
  (J)  jump      : |nu| - max(|lam|,|mu|) <= 1          [Q484, measured, fitted]
  (U)  uncovered : |nu \ (lam u mu)| <= 1               [mechanism, derived]
Enumerate triples permitted by the support bound nu>=lam, nu>=mu (Motegi's order).
"""
import sys; sys.path.insert(0,'/home/clio/projects/scratch')
import q484_engine as E
from itertools import combinations

def crits(N):
    S=sorted(E.subsets(N), key=lambda M:(len(M),M))
    tr=[]
    for lam in S:
        for mu in S:
            for nu in S:
                if E.geq(nu,lam) and E.geq(nu,mu): tr.append((lam,mu,nu))
    return tr

for N in (2,3,4,5):
    tr=crits(N)
    dis=[]
    for lam,mu,nu in tr:
        J = len(nu)-max(len(lam),len(mu)) <= 1
        U = len(set(nu)-(set(lam)|set(mu))) <= 1
        if J!=U: dis.append((lam,mu,nu,J,U))
    print(f'N={N}: {len(tr)} triples permitted by support bound; disagreements J vs U: {len(dis)}')
    for d in dis[:12]:
        lam,mu,nu,J,U=d
        print(f'   lam={lam} mu={mu} nu={nu}   jump-crit Laurent={J}   uncovered-crit Laurent={U}')
    if len(dis)>12: print(f'   ... and {len(dis)-12} more')
