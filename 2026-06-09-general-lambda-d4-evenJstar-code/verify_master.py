"""
Cross-check: recompute G_lambda(i) two ways and confirm they agree.
  (A) M_j expansion:  G = sum_j C(m,j) i^j (1+i)^j M_j           [job1 G_gaussian]
  (B) master identity: G = <s_lambda, psi^m>,  psi = h1^2 + i(1+i) e2,
      computed in the power-sum basis over Z[i] (symfunc.py).
"""
from sympy import I, Integer, expand
import symfunc as sf
from job1_tie_census import partitions, M_vector, G_gaussian

def psi():
    # psi = h1^2 + i(1+i) e2  ; h1=p1, e2 in p-basis. coeff i(1+i) = i-1.
    h1sq = sf.power(sf.h1(), 2)
    e2 = sf.e_n(2)
    coeff = I*(1+I)            # = -1 + i
    term2 = {k: v*coeff for k,v in e2.items()}
    return sf.add(h1sq, term2)

def G_master(lam, m):
    ps = psi()
    pw = sf.power(ps, m)
    val = sf.inner(sf.schur_p(lam), pw)
    val = expand(val)
    re = val.as_real_imag()[0]; im = val.as_real_imag()[1]
    return int(re), int(im)

if __name__=="__main__":
    nbad=0; ntot=0
    for m in range(1,8):
        for lam in partitions(2*m):
            Ms=M_vector(lam,m)
            reA,imA=G_gaussian(lam,m,Ms)
            reB,imB=G_master(lam,m)
            ntot+=1
            if (reA,imA)!=(reB,imB):
                nbad+=1
                print("MISMATCH",lam,m,(reA,imA),(reB,imB))
    print(f"verified {ntot} shapes (m<=7): {ntot-nbad} agree, {nbad} mismatch")
    # also confirm (2,2) is exactly 0 both ways
    print("(2,2):", G_gaussian((2,2),2,M_vector((2,2),2)), G_master((2,2),2))
