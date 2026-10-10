import sympy as sp
from core import Engine
from warn import G, t
from walks import AB_from_walks

def report(k, ell, n):
    w = 2*ell+2; E = Engine(n)
    Gt = sp.expand(G(E, k, ell, t) * E.prodx**k)
    P = sp.Poly(Gt, *E.x)
    out = dict(k=k,ell=ell,w=w,n=n,nmono=0,strict=[],anysign=[])
    for mono, coeff in zip(P.monoms(), P.coeffs()):
        co = sp.expand(coeff)
        d = {int(kk[0]):int(v) for kk,v in (sp.Poly(co,t).as_dict().items() if co.has(t) else {(0,):co}.items())}
        l1  = sum(abs(v) for v in d.values())
        pos = sum(v for v in d.values() if v>0); neg = sum(-v for v in d.values() if v<0)
        A,B,TOT = AB_from_walks(k,w,mono)
        out['nmono'] += 1
        assert A-B == int(co.subs(t,1))
        if pos>A or neg>B: out['strict'].append((max(pos-A,neg-B),mono,d,A,B,TOT))
        if l1>TOT:          out['anysign'].append((l1-TOT,mono,d,A,B,TOT))
    out['strict'].sort(reverse=True); out['anysign'].sort(reverse=True)
    return out
