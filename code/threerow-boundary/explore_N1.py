import sympy as sp
b,P = sp.symbols('b P', integer=True)
N1s = 2*P*b**2 + 10*P*b + 12*P + 7*b**2 + 17*b + 12

# express around (P+2)
for form,expr in [
  ("2(P+2)(b+2)(b+3)+3(b^2-b-4)", 2*(P+2)*(b+2)*(b+3) + 3*(b**2-b-4)),
  ("2(P+1)(b+2)(b+3)+ ?", None),
]:
    if expr is None: continue
    print(form, "diff=", sp.simplify(sp.expand(expr)-N1s))

# reduce N1 mod (P+1) and (P+2): substitute P=-1, P=-2
print("N1 at P=-1 (mod P+1):", sp.factor(N1s.subs(P,-1)))
print("N1 at P=-2 (mod P+2):", sp.factor(N1s.subs(P,-2)))
# so N1 = (P+1)*Q1 + r1 ?
Q1, r1 = sp.div(sp.Poly(N1s,P), sp.Poly(P+1,P))
print("N1 = (P+1)*(", sp.factor(Q1.as_expr()), ") +", sp.factor(r1.as_expr()))
Q2, r2 = sp.div(sp.Poly(N1s,P), sp.Poly(P+2,P))
print("N1 = (P+2)*(", sp.factor(Q2.as_expr()), ") +", sp.factor(r2.as_expr()))
