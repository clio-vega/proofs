from math import comb
def v2(n):
    n=int(n)
    if n==0: return 10**9
    r=0
    while n%2==0: n//=2; r+=1
    return r
def C(n,k):
    return comb(n,k) if 0<=k<=n else 0
def falling(j,k):
    p=1
    for i in range(k): p*=(j-i)
    return p
def H(a,b,j):
    return (a+3)*(a+4)*(b+2)*(b+3) -6*(a+3)*(b+2)*C(j,1) -6*(a*b+a+2*b)*C(j,2) +36*C(j,3)+72*C(j,4)
def Q3(a,b,j):
    return (a-1)*(b-2)*H(a,b,j) - 720*C(j,6)

# Test split: T1=(a-1)(b-2)C(b+3,j)H, T2 = C(b+3,j)*falling(j,6) ; check both >= v2(b+3)+1 (a even)
bad_full=0; bad_T1=0; bad_T2=0; n=0
for a in range(2,60,2):      # a even
    for b in range(2,a+1):    # b<=a typical; partition rows a>=b
        for j in range(6, b+4):  # tip active region j>=6, j<=b+3
            cb=C(b+3,j)
            if cb==0: continue
            n+=1
            full=v2(cb)+v2(Q3(a,b,j))
            target=v2(b+3)+1
            T1=v2(a-1)+v2(b-2)+v2(cb)+v2(H(a,b,j))
            T2=v2(cb)+v2(falling(j,6))
            if full<target: bad_full+=1
            if T1<target: bad_T1+=1
            if T2<target: bad_T2+=1
print(f"tested {n} (a even, j>=6). full<target:{bad_full}  T1<target:{bad_T1}  T2<target:{bad_T2}")

# also j<6 region (tip=0): Q3=(a-1)(b-2)H
bad_low=0; nl=0
for a in range(2,60,2):
    for b in range(2,a+1):
        for j in range(0,6):
            cb=C(b+3,j)
            if cb==0: continue
            nl+=1
            if v2(cb)+v2(Q3(a,b,j))<v2(b+3)+1: bad_low+=1
print(f"j<6 region tested {nl}, failures {bad_low}")
