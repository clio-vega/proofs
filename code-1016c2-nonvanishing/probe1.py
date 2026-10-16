import sys; sys.path.insert(0,'.')
import sympy as sp
import q484_engine as E
q=E.q
from gate_q0 import build, solve_c, linext, qval, limit_q0
w,S,A = build(2,'C')
c = solve_c(S,A,(),())
f = c[(1,2)]
print('c^{12}_{00} =', sp.factor(sp.cancel(f)))
num,den = sp.fraction(sp.cancel(f))
print('\nnum =', sp.factor(num))
print('\nden =', sp.factor(den))
print('\nqval(c) =', qval(f))
s=sp.Symbol('s')
cs = sp.cancel(f.subs(w[1], s/w[0]))
print('\nc with w2=s/w1 :', sp.factor(cs))
lead = sp.cancel(sp.together((s-1)*cs))
print('\n(s-1)*c =', sp.factor(lead))
n1,d1 = sp.fraction(lead)
print('\n  at s=1: num=',sp.factor(sp.expand(n1).subs(s,1)),' den=',sp.factor(sp.expand(d1).subs(s,1)))
res = sp.cancel(sp.expand(n1).subs(s,1)/sp.expand(d1).subs(s,1))
print('\nRESIDUE =', sp.factor(res), '   qval =', qval(res))
print('\nc at q=0 :', sp.factor(limit_q0(f)[0]), '  (qval', limit_q0(f)[1],')')
