import sys; sys.path.insert(0,'.')
from fibre import enumerate_fibre
import sympy as sp, re
t=sp.Symbol('t')
rows=[]
for line in open('ground.txt'):
    m=re.match(r'k=(\d+) ell=(\d+) h=(\d+) w=(\d+) n=(\d+): c=\((.*)\)  c\(1\)=(\S+)', line)
    if m:
        k,ell,h,w,n=map(int,m.groups()[:5]); c=sp.expand(sp.sympify(m.group(6)))
        rows.append((k,ell,h,w,n,c))
bad=0
print(" k ell  n |  |fibre|  deg_t  ord_t=1 |  Warnaar c(1)   sum c^-   match")
for (k,ell,h,w,n,c) in rows:
    F=enumerate_fibre(h,w,n); s=sum(f[2] for f in F)
    c1=sp.simplify(c.subs(t,1))
    if c==0: deg,o='-','-'
    else:
        deg=sp.Poly(c,t).degree(); o=0
        while sp.simplify(c.subs(t,1)) == 0 and o<30:
            c2=sp.cancel(c/(t-1)**(o+1))
            if sp.simplify(c2.subs(t,1))!=0: o+=1; break
            o+=1
        o = sp.Poly(c,t).all_coeffs() and sum(1 for _ in iter(int,1)) if False else o
    ok = (sp.simplify(c1-s)==0)
    if not ok: bad+=1
    print(f"{k:>2}{ell:>4}{n:>3} | {len(F):>8} {str(deg):>6} {str(o):>8} | {str(c1):>13} {s:>9}   {ok}")
print("MISMATCHES:", bad, " / rows:", len(rows))
