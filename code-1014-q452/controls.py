import sys; sys.path.insert(0,'.')
from fibre import gvars, sigma_symbolic, enumerate_fibre
import sympy as sp
t=sp.Symbol('t')

print("CONTROL 1 (Q430 reproduction): index-aware weight hitting 1-t on the ell=0, alpha=(1,1) fibre")
for k in [1,2,3,4]:
    h,w=2*k,2; g=gvars(h,w)
    sub={g[(h,2)]:1, g[(1,1)]:1, g[(h,1)]:t}
    S=sigma_symbolic(h,w,2,g).subs(sub)
    F=enumerate_fibre(h,w,2)
    print(f"   k={k} h={h}: |fibre|={len(F)}  Sigma_2 = {sp.expand(S)}   == 1-t ? {sp.expand(S-(1-t))==0}")

print()
print("CONTROL 2 (satisfiability at (k,ell)=(1,0)): a POLYNOMIAL index-aware weight matching c_2 AND c_4")
h,w=2,2; g=gvars(h,w)
sol={g[(1,1)]:2-t, g[(1,2)]:-t**2, g[(2,1)]:sp.Integer(1), g[(2,2)]:sp.Integer(1)}
tg={2: 1-t, 4: sp.expand((t-1)**2*(t+2)), 6: sp.expand(-(t-1)**3*(t**3+3*t**2+6*t+5))}
print(f"   weight: g_1_1={2-t}, g_1_2={-t**2}, g_2_1=1, g_2_2=1   (all in QQ[t])")
for n in (2,4,6):
    S=sp.expand(sigma_symbolic(h,w,n,g).subs(sol))
    print(f"   n={n}: Sigma={sp.factor(S)}   target={sp.factor(tg[n])}   MATCH? {sp.expand(S-tg[n])==0}"
          + ("" if sp.expand(S-tg[n])==0 else f"   residual={sp.factor(sp.expand(tg[n]-S))}"))

print()
print("CONTROL 3 (satisfiability IN DOMAIN, (k,ell)=(2,0)): solve n=2,4,6 at h=4,w=2")
h,w=4,2; g=gvars(h,w)
tg2={2: 1-t, 4: sp.expand(3*(t-1)**2), 6: sp.expand((t-1)**3*(t**2-2*t-14))}
# ansatz: kill the coupling by setting several g to simple values, solve the rest
a,b,c = sp.symbols('a b c')
ans={g[(4,2)]:sp.Integer(1), g[(4,1)]:sp.Integer(0), g[(1,2)]:sp.Integer(0),
     g[(2,2)]:sp.Integer(0), g[(3,2)]:sp.Integer(0),
     g[(1,1)]:a, g[(2,1)]:b, g[(3,1)]:c}
eqs=[sp.expand(sigma_symbolic(h,w,n,g).subs(ans)-tg2[n]) for n in (2,4,6)]
for n,e in zip((2,4,6),eqs): print(f"   n={n}: {sp.factor(e)} = 0")
sols=sp.solve(eqs,[a,b,c],dict=True)
print("   solutions:", sols)
for s in sols:
    full={**{kk:vv.subs(s) if hasattr(vv,'subs') else vv for kk,vv in ans.items()}}
    ok=[sp.expand(sigma_symbolic(h,w,n,g).subs(full)-tg2[n])==0 for n in (2,4,6)]
    print(f"   verify n=2,4,6: {ok}")
