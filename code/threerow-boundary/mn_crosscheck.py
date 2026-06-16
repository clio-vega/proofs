from mn import Mj
# verify M_{b+1},M_{b+2},M_{b+3} closed forms against Murnaghan-Nakayama
bad=[]
for m in range(6,30):
    for b in range(3,m):
        a=2*m-3-b
        if a<b: continue
        c=3
        N2=a*(b+3)-(b**2+2*b+3)
        N1=a*b**2+5*a*b+6*a-b**3-b**2-4*b-6
        Mb3=(b-2)*(b+2)*(b+3)//6
        Mb2=(b-2)*(b+2)*N2//6
        Mb1=(b-2)*(a-b+1)*N1//12
        for i,pred in [(1,Mb1),(2,Mb2),(3,Mb3)]:
            j=b+i
            if j>m: continue
            actual=Mj((a,b,c),j,m)
            if actual!=pred: bad.append(((a,b,c),m,i,actual,pred))
print("M closed-form mismatches:", bad[:10], "count", len(bad))
