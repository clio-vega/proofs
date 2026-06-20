import sympy as sp, pickle
a,b,j=sp.symbols('a b j')
d=pickle.load(open('/home/clio/projects/code/threerow-c5/H5sym.pkl','rb'))
hcoeffs=[sp.sympify(s) for s in d['hcoeffs']]
def binom_poly(jj,k):
    r=sp.Integer(1)
    for t in range(k): r*=(jj-t)
    return r/sp.factorial(k)
H5=lambda J: sp.expand(sum(hcoeffs[k]*binom_poly(J,k) for k in range(9)))
for J in [0,1,2,3,4,5]:
    print(f"H5({J}) =", sp.factor(H5(J)))
