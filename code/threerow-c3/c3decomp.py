import sympy as sp
a,b,j = sp.symbols('a b j', integer=True)
K3=(a-1)*(a+3)*(a+4)*(b-2)*(b+2)*(b+3)
pert=(-j**6+15*j**5 + j**4*(3*a*b-6*a-3*b-79) + j**3*(-12*a*b+24*a+12*b+201)
      + j**2*(-3*a**2*b**2+3*a**2*b+6*a**2-3*a*b**2+24*a*b-36*a+6*b**2-27*b-244)
      + j*(-3*a**2*b**2-3*a**2*b+18*a**2-9*a*b**2-15*a*b+66*a+12*b**2+18*b+36))

# express pert in basis C(j,k), k=1..6 : pert = sum_k coef_k * binomial(j,k)
# coef_k = k-th forward difference at 0 = sum_{i=0}^k (-1)^{k-i} C(k,i) pert(j=i)
def fdiff(k):
    s=sp.Integer(0)
    for i in range(k+1):
        s+= (-1)**(k-i)*sp.binomial(k,i)*pert.subs(j,i)
    return sp.expand(s)
print("pert = sum_k coef_k * C(j,k):")
for k in range(1,7):
    print(f"  coef_{k} =", sp.factor(fdiff(k)))

# Also K3 parity-relevant factorization already known.
# Now: the 'inhomogeneous tip' P6 candidate. coef_6 should be 6!*(-1)= -720 => -720 C(j,6) leading
print("\nNow Q3 itself in binomial basis (Q3 = K3 + pert), constant = K3:")
Q3=K3+pert
def fdiffQ(k):
    s=sp.Integer(0)
    for i in range(k+1):
        s+= (-1)**(k-i)*sp.binomial(k,i)*Q3.subs(j,i)
    return sp.expand(s)
for k in range(0,7):
    print(f"  C(j,{k}) coef =", sp.factor(fdiffQ(k)))
