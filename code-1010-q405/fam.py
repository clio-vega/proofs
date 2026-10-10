import sympy as sp
from core import Engine
from warn import G, t
from walks import AB_from_walks
out=[]
for ell in range(0,5):
    M=ell+2; n=2*M; w=2*ell+2; E=Engine(n)
    Gt=sp.expand(G(E,1,ell,t)*E.prodx)
    c=sp.expand(sp.Poly(Gt,*E.x).as_dict().get(tuple([1]*n),0))
    CM=int(sp.binomial(2*M,M)-sp.binomial(2*M,M-1))
    pred=CM-(2*M-1)*t**(M-1)+t**(M+1)
    A,B,TOT=AB_from_walks(1,w,tuple([1]*n))
    L1=sum(abs(int(v)) for v in sp.Poly(c,t).as_dict().values())
    print('ell=%d M=%d n=N=%d  c=%-30s match=%s | A=%d(C_M=%d) B=%d(2M-2=%d) TOT=%d L1=%d | (S)fails=%s (S+-)fails=%s'
          %(ell,M,n,str(c),sp.expand(c-pred)==0,A,CM,B,2*M-2,TOT,L1,(2*M-1)>B,L1>TOT), flush=True)
