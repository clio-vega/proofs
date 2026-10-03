"""Cases 2 and 3 of Theorem M returned IDENTICAL counts. Conjecture: the involution
   y |-> -y  composed with swapping the roles of y and y'
is a bijection of the test set that exchanges Case 2 and Case 3.
On the hyperplane sum y = sigma:  k(-y) = sum_i (y_i)_+ = k(y) + sigma,
so {y in B : sum=sigma, k<=j}  ->  {z in -B : sum=-sigma, k<=j+sigma}  bijectively.
Verify both halves explicitly."""
import itertools
def k(y): return sum(max(-t,0) for t in y)
ok=bad=0; casemap={}
for m,lo,hi in [(2,-3,3),(3,-2,2)]:
    rngs=[(P,Q) for P in range(lo,hi+1) for Q in range(P,hi+1)]
    for B in itertools.product(rngs,repeat=m):
        nB=tuple((-Q,-P) for (P,Q) in B)
        for sigma in range(sum(P for P,_ in B), sum(Q for _,Q in B)+1):
            Y=[y for y in itertools.product(*[range(P,Q+1) for (P,Q) in B]) if sum(y)==sigma]
            if len(Y)<2: continue
            for j in sorted({k(y) for y in Y}):
                S=[y for y in Y if k(y)<=j]
                if len(S)<2: continue
                # half 1: the image is exactly the sublevel set at level j+sigma over -B
                Z=[z for z in itertools.product(*[range(P,Q+1) for (P,Q) in nB]) if sum(z)==-sigma]
                Simg=sorted(tuple(-t for t in y) for y in S)
                Starget=sorted(z for z in Z if k(z)<=j+sigma)
                if Simg==Starget: ok+=1
                else: bad+=1
                # half 2: case labels are exchanged
                for y in S:
                    for yp in S:
                        for i in range(m):
                            if y[i]<=yp[i]: continue
                            c = 1 if (y[i]>=1 and yp[i]<=-1) else (2 if y[i]>=1 else 3)
                            z, zp = tuple(-t for t in yp), tuple(-t for t in y)   # roles swapped
                            assert z[i]>zp[i]
                            c2 = 1 if (z[i]>=1 and zp[i]<=-1) else (2 if z[i]>=1 else 3)
                            casemap[(c,c2)]=casemap.get((c,c2),0)+1
print(f" sublevel-set image matches target: {ok} ok, {bad} mismatches")
print(" case -> image case:",dict(sorted(casemap.items())))
