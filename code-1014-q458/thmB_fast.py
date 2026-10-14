import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1014-q458')
from censusk import c_1n
ok=fail=0; rows=[]
for k in range(1,5):
    for ell in range(0,4):
        M=k+ell+1; N=2*M
        nmax={1:36,2:30,3:20,4:14}[k]
        for n in range(2,nmax+1,2):
            if n%N and (n+2*k)%N: continue
            sign = +1 if n%N==0 else -1
            D = n*(n+2*k)//(4*M)
            d=c_1n(n,k,ell)
            gd,gc=max(d),d[max(d)]
            good=(gd==D and gc==sign); ok+=good; fail+=not good
            rows.append((k,ell,n,'B' if sign>0 else 'A',D,sign,gd,gc,good))
            if not good: print("MISMATCH",rows[-1],flush=True)
print(f"THEOREM B (k<=4): ok={ok} fail={fail}")
print("sample (k,ell,n,branch,pred_deg,pred_coef,got_deg,got_coef,ok):")
for r in rows[::max(1,len(rows)//20)]: print("  ",r)
