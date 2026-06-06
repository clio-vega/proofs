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
# Test: does mu_lb_{k+1} = (2 a^2 + b^2/mu_k)/(k+1)^2  >  R_{k+1}(x_k) ?
# i.e. does the S>=0 lower bound on mu already force x_{k+1}<1, given TRUE mu_k,x_k?
bad=0; worst=mp.mpf('1e9'); winfo=None
for m in range(4,160):
    N,Q=seqs(m)
    mu=[None]+[N[k]/N[k-1] for k in range(1,m+1)]
    x=[None]+[ (2*m-k+1)*Q[k]/((m-k)*N[k]) if k<m else None for k in range(1,m+1)]
    for k in range(1,m-1):  # predict x_{k+1}, need k+1<=m-1 => k<=m-2
        a=m-k; b=2*m-k+1
        mulb=(2*a*a + b*b/mu[k])/(k+1)**2
        a1=m-k-1; b1=2*m-k
        R = b1*a*(1-x[k])/(a1*(k+1))
        margin = mulb - R   # want >0
        if margin< -mp.mpf('1e-30'): bad+=1
        if margin<worst: worst=margin; winfo=(m,k,m-k,float(mulb),float(R))
print("cases where S>=0 lower bound FAILS to give x_{k+1}<1:", bad)
print("worst margin (mulb-R):", float(worst)," at (m,k,j=m-k,mulb,R)=",winfo)
