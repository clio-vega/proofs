"""Heavy piece (star1):  v2[ C(a+4,j) C(b+3,j) H(j) ] >= v2H(0) + s2(j) - (j+3)/2,
   a odd, b even, 4<=j<=b.  v2H(0)=v2(a+3)+v2(b+2).
Find H's heavy/tip decomposition and probe sub-bounds.
"""
import sympy as sp
a,b,j=sp.symbols('a b j')
H=( (a+3)*(a+4)*(b+2)*(b+3) -6*(a+3)*(b+2)*sp.binomial(j,1)
    -6*(a*b+a+2*b)*sp.binomial(j,2)+36*sp.binomial(j,3)+72*sp.binomial(j,4) )
# Try H = (a+3)(b+2)*G - tip ?  Compute H - (a+3)(b+2)*[(a+4)(b+3)-6 C(j,1)]
rem = sp.expand(H - (a+3)*(b+2)*((a+4)*(b+3)-6*sp.binomial(j,1)))
print("H - (a+3)(b+2)[(a+4)(b+3)-6C(j,1)] =")
print(sp.factor(rem))
print("  in binomial basis:", sp.expand(rem))
# Express rem in terms of C(j,2),C(j,3),C(j,4)
c2,c3,c4=sp.binomial(j,2),sp.binomial(j,3),sp.binomial(j,4)
print("\nrem coeffs: -6(ab+a+2b)C(j,2)+36C(j,3)+72C(j,4):")
print("  ab+a+2b =", sp.factor(a*b+a+2*b))
# alt: ab+a+2b = (a+2)(b+1)-2
print("  (a+2)(b+1)-2 - (ab+a+2b) =", sp.expand((a+2)*(b+1)-2-(a*b+a+2*b)))
