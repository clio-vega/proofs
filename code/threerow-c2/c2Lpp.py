from math import comb
def v2(n):
    n=abs(int(n))
    if n==0:return 10**9
    k=0
    while n%2==0:n//=2;k+=1
    return k
def vC(n,k):
    if k<0 or k>n:return 10**9
    return v2(comb(n,k))
def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j)+j*(j-1)*(j-2)*(j-3)

# (L''): v2C(O,j)+v2Q(j) >= v2F+1, O=odd one of {a+3,b+2}, F=even one of {a,b-1}
# a,b even: a+3 odd=O, b+2 even=E; a even=F, b-1 odd. so O=a+3, F=a.
# a,b odd:  a+3 even=E, b+2 odd=O; a odd, b-1 even=F. so O=b+2, F=b-1.
badLpp=0; ex=[]
# also test full (L) again as control
badL=0
for A in range(3,400):
    for B in range(2,A+1):
        if (A+B)%2: continue
        if A%2==0:  # both even
            O=A+3; F=A
        else:       # both odd
            O=B+2; F=B-1
        for j in range(1,B+1):
            lhs=vC(O,j)+v2(Q(A,B,j)); rhs=v2(F)+1
            if lhs<rhs:
                badLpp+=1
                if len(ex)<12: ex.append((A,B,j,lhs,rhs))
            # full L:
            T=vC(A+3,j)+vC(B+2,j)+v2(Q(A,B,j))-v2(Q(A,B,0))
            if T<1-v2(j): badL+=1
print("(L'') v2C(O,j)+v2Q(j) >= v2F+1  FAILURES:",badLpp)
print("  examples:",ex[:12])
print("(L) control failures:",badL)
