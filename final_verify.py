import mpmath as mp
mp.mp.dps=40
def seqs(m):
    poly=[mp.mpc(1),mp.mpc(1,1),mp.mpc(1)]; res=[mp.mpc(1)]
    for _ in range(m):
        new=[mp.mpc(0)]*(len(res)+2)
        for a,ca in enumerate(res):
            for b,cb in enumerate(poly): new[a+b]+=ca*cb
        res=new
    r=res; N=[abs(x)**2 for x in r]
    Q=[None]+[(r[k]*mp.conj(r[k-1])).imag for k in range(1,len(r))]
    return N,Q
# check limiting recursion y_{j-1}=j/(j-1)(1-y_j) against exact x near top, large m
m=300; N,Q=seqs(m)
def x(k): a=m-k;b=2*m-k+1; return float(b*Q[k]/(a*N[k]))
print("near-top: j, exact y_j=x_{m-j}, predicted j/(j-1)*(1-y_j) vs actual y_{j-1}")
for j in range(8,1,-1):
    yj=x(m-j); yjm1=x(m-(j-1))
    pred=j/(j-1)*(1-yj)
    print(f" j={j}: y_j={yj:.6f}  pred y_{{j-1}}={pred:.6f}  actual={yjm1:.6f}  err={pred-yjm1:+.2e}")
# confirm invariant band 1/2<y_j<(j+1)/2j holds for these large-m top values
print("\n band check 1/2 < y_j < (j+1)/(2j) for j=1..8 at m=300:")
for j in range(1,9):
    yj=x(m-j); print(f" j={j}: {0.5:.4f} < {yj:.6f} < {(j+1)/(2*j):.6f} : {0.5<yj<(j+1)/(2*j)}")
