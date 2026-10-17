"""Is there a SHORT-root collision?  For m in M and eps^2=1, is
     Atilde_{M,I} |_{w_m=eps}  =  sigma * Atilde_{M\{m},I} |_{w_m=eps}
with sigma independent of M and I?  N=3, type C, symbolic q.
"""
import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1017c2-coboundary')
import sympy as sp
import q484_engine as E
q=E.q
N=3; w=E.wvars(N); S=sorted(E.subsets(N),key=lambda M:(len(M),M))
A={(M,I):E.Atilde(N,w,M,I,'C') for M in S for I in S}
for eps in (1,-1):
    print(f'===== w_m = {eps} =====',flush=True)
    for m in range(1,N+1):
        vals={}
        for M in S:
            if m not in M: continue
            Mm=tuple(j for j in M if j!=m)
            for I in S:
                top=sp.cancel(A[(M,I)]); bot=sp.cancel(A[(Mm,I)])
                t=sp.cancel(top.subs(w[m-1],eps)); b=sp.cancel(bot.subs(w[m-1],eps))
                if b==0:
                    tag='top=0,bot=0' if t==0 else 'BOT=0,TOP!=0'
                    vals.setdefault(tag,[]).append((M,I)); continue
                r=sp.cancel(sp.simplify(t/b))
                vals.setdefault(sp.sstr(sp.factor(r)),[]).append((M,I))
        print(f'  m={m}: {len(vals)} distinct ratio(s)',flush=True)
        for k,v in sorted(vals.items(), key=lambda kv:-len(kv[1])):
            print(f'     {len(v):3d}x  {k}',flush=True)
