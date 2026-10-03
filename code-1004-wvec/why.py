"""Why did all three stratifications give IDENTICAL counts (11580)?
On a slice sum_i y_i = sigma is CONSTANT, so
    sum_i (-y_i)_+ = (||y||_1 - sigma)/2,   sum_i (y_i)_+ = (||y||_1 + sigma)/2.
All three are affine functions of each other => they induce the SAME stratification.
My 'refusal tests' were the same predicate three times.  Verify that, then build real ones."""
import gen, tent
def ys(mu,lam,n,m,b):
    t=sum(mu)+b; out=[]
    for nu in gen.box(mu,n,m):
        if sum(nu)!=t: continue
        lr=gen.LR(nu,lam,n,m)
        if any(R-L+1<=0 for (L,R) in lr): continue
        out.append(tuple(nu[i]-tent.lam_at(lam,i-1,n,m)-1 for i in range(m)))
    return out
same=diff=0
for m in (2,3,4):
  for n in range(m,m+5):
    for (mu,lam) in gen.pairs(n,m,7):
        d=sum(lam)-sum(mu)
        for b in range(0,d+1):
            Y=ys(mu,lam,n,m,b)
            if len(Y)<2: continue
            # do the three functions induce the same ordered partition of Y?
            def strat(f):
                v=sorted({f(y) for y in Y})
                return tuple(frozenset(y for y in Y if f(y)<=j) for j in v)
            a=strat(lambda y: sum(max(-t,0) for t in y))
            b2=strat(lambda y: sum(max(t,0) for t in y))
            c=strat(lambda y: sum(abs(t) for t in y))
            if a==b2==c: same+=1
            else: diff+=1
print(f" slices where the three stratifications coincide: {same}; differ: {diff}")
