import sympy as sp
from ribbon import *
from sympy import Matrix, eye, simplify, symbols
t,z = sp.Symbol('t'), sp.Symbol('z')

print("=== GATE 1a: is Z_{lam,mu} = <s_lam, R_e(t) s_mu> square? ===")
for e in (2,3):
    for N in (6,):
        rows=set(); cols=set()
        for mu in all_parts_upto(N-e):
            for lam,ht in add_ribbons(mu,e,N+3):
                if sum(lam)<=N: rows.add(lam); cols.add(mu)
        print(f"  e={e}, |lam|,|mu| <= {N}: support has |lam|=|mu|+{e} always:",
              all(sum(l)==sum(m)+e for m in all_parts_upto(N-e) for l,_ in add_ribbons(m,e,N+3)))
        print(f"    #rows hit={len(rows)}  #cols hit={len(cols)}  -> rectangular, Z^-1 undefined")

print()
print("=== GATE 1b: w_e(lam') == w_e(lam)?  and core_e(lam')==core_e(lam)' ? ===")
bad=0;tot=0
for e in (2,3,4):
    for n in range(0,11):
        for lam in partitions(n):
            c1,w1=ecore_weight(lam,e); c2,w2=ecore_weight(conj(lam),e)
            tot+=1
            if w1!=w2 or c2!=conj(c1): bad+=1; print("  FAIL",e,lam)
print(f"  {tot-bad}/{tot} partitions: w_e and core commute with conjugation")
