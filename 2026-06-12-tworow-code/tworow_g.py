from math import comb as C

def v2(n):
    if n==0: return float('inf')
    v=0
    while n%2==0:
        n//=2; v+=1
    return v
def s2(n):
    return bin(n).count('1')

def Mj(m,a,b,j):
    return C(2*(m-j), b-j) - (C(2*(m-j), b-j-1) if b-j-1>=0 else 0)

ok_g=True
ok_tie=True
ok_subset=True
ties=[]
for m in range(1,15):
    for b in range(0,m+1):
        a=2*m-b
        # compute val(j) directly
        vals={}
        for j in range(0,b+1):
            Mjv=Mj(m,a,b,j)
            cm = C(m,j)*Mjv
            if cm==0:
                continue
            vals[j]=j+2*v2(cm)
        V=min(vals.values())
        Jstar=sorted([j for j in vals if vals[j]==V])
        # g(j) formula check
        for j in vals:
            g_direct = vals[j]-vals[0]
            g_formula = j + 2*(v2(C(b,j)) + v2(C(a+1,j)) - s2(j))
            if g_direct != g_formula:
                ok_g=False; print("g MISMATCH",m,a,b,j,g_direct,g_formula)
        # subset check
        if not set(Jstar).issubset({0,2}):
            ok_subset=False; print("SUBSET FAIL",m,a,b,Jstar)
        if 0 not in Jstar:
            ok_subset=False; print("0 not in Jstar",m,a,b,Jstar)
        # tie characterization: Jstar={0,2} iff C(b,2),C(a+1,2) both odd
        is_tie = (Jstar==[0,2])
        pred_tie = (b>=2) and (C(b,2)%2==1) and (C(a+1,2)%2==1)
        if is_tie!=pred_tie:
            ok_tie=False; print("TIE MISMATCH",m,a,b,Jstar,"pred",pred_tie)
        if is_tie:
            ties.append((m,a,b))
print("g(j) formula correct (m<=14):", ok_g)
print("J* subset {0,2}, 0 in J* (m<=14):", ok_subset)
print("tie char C(b,2),C(a+1,2) both odd (m<=14):", ok_tie)
print("number of ties:", len(ties))
# also verify mod-4 description: tie iff b in{2,3 mod4} and a in {1,2 mod4}
ok_mod=True
for m in range(1,15):
    for b in range(2,m+1):
        a=2*m-b
        Jstar_tie = (C(b,2)%2==1 and C(a+1,2)%2==1)
        pred = (b%4 in (2,3)) and (a%4 in (1,2))
        if Jstar_tie!=pred:
            ok_mod=False; print("MOD MISMATCH",m,a,b)
print("mod-4 tie description correct:", ok_mod)
print("sample ties (m,a,b):", ties[:12])
