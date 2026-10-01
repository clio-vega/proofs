"""Verify the structural identities and the radial/tail decomposition."""
from regions import *
from fractions import Fraction as F

def data(n,A,B,u,Tm,Tp):
    C,D = B+1, A+n-1
    sigma0 = F(u+B+D,2); t0 = F(u+B-D,2)
    rows=[]
    for t in range(Tm,Tp+1):
        at,bt,ct,dt = endpoints(t,n,A,B,u)
        nt, mt = bt-at+1, dt-ct+1
        rows.append(dict(t=t,at=at,bt=bt,ct=ct,dt=dt,nt=nt,mt=mt,
                         L=F(nt+mt-2,2), delta=F(abs(nt-mt),2)))
    return C,D,sigma0,t0,rows

def check_identities(n,A,B,u,Tm,Tp):
    C,D,s0,t0,rows = data(n,A,B,u,Tm,Tp)
    for R in rows:
        t=R['t']
        if R['at']+R['dt'] != t+D: return ('ad',t)
        if R['bt']+R['ct'] != u+B-t: return ('bc',t)
        if R['nt']<1 or R['mt']<1: continue
        if F(R['at']+R['ct']+R['bt']+R['dt'],2) != s0: return ('sigma',t)
        if R['delta'] != abs(F(t)-t0): return ('delta',t)
    # L concave & 1-Lipschitz in t over ALL t (formula-level, all t in Z)
    Lf = lambda t: F(sum([min(u-1-t,B)-max(t,A)+1, min(t+n-1,D)-max(u-t,C)+1])-2,2)
    for t in range(-30,40):
        if 2*Lf(t) < Lf(t-1)+Lf(t+1): return ('Lconc',t)
        if abs(Lf(t+1)-Lf(t))>1: return ('Llip',t)
    return None

def radial(n,A,B,u,Tm,Tp):
    """Return (sigma0, beta dict r->count, gamma dict y->tail, G from gamma)."""
    C,D,s0,t0,rows = data(n,A,B,u,Tm,Tp)
    beta={}
    for R in rows:
        if R['nt']<1 or R['mt']<1: continue
        r = R['delta']
        while r <= R['L']:
            beta[r]=beta.get(r,0)+1
            r = r+1
    if not beta: return s0,beta,{},{}
    rmax=max(beta); rmin=min(beta)
    gamma={}; acc=0
    r=rmax
    # gamma must be built down to the BOTTOM of the lattice (r>=0), not to min(beta):
    # the lattice of |s-sigma0| starts at 0 or 1/2, which may be below min(beta).
    bot = rmin - (rmin if rmin==int(rmin) else rmin)   # keep lattice class
    stop = rmin
    while stop - 1 >= 0: stop = stop - 1
    while r>=stop:
        acc+=beta.get(r,0); gamma[r]=acc; r=r-1
    Gr={}
    for y,v in gamma.items():
        if v>0:
            Gr[s0+y]=v; Gr[s0-y]=v
    return s0,beta,gamma,{int(k):v for k,v in Gr.items()}

def beta_is_truncated_concave(beta):
    """beta nonneg, interval support, and increments non-increasing ON support."""
    if not beta: return True
    ks=sorted(beta); lo,hi=ks[0],ks[-1]
    seq=[]; r=lo
    while r<=hi: seq.append(beta.get(r,0)); r=r+1
    nz=[i for i,v in enumerate(seq) if v>0]
    if nz and nz[-1]-nz[0]+1!=len(nz): return False
    sup=seq[nz[0]:nz[-1]+1] if nz else []
    inc=[sup[i+1]-sup[i] for i in range(len(sup)-1)]
    return all(inc[i]>=inc[i+1] for i in range(len(inc)-1))
