import sys; sys.path.insert(0,'/home/clio/projects/code/threerow-c1')
from mn import Mj
from math import comb
from fractions import Fraction

def Q(a,b,j): return a*(b-1)*((a+3)*(b+2)-2*j*j) + j*(j-1)*(j-2)*(j-3)
def Mj_formula(a,b,j):
    m=(a+b+2)//2; N=2*(m-j); k=b-j
    if k<0: return None
    cb=comb(N,k)
    return Fraction(cb*(a-b+1)*Q(a,b,j), 2*(a+3-j)*(b+2-j)*(b+1-j))

bad=0; tested=0
for A in range(3,24):
    for B in range(2,A+1):   # b>=c=2
        if (A+B)%2!=0: continue
        m=(A+B+2)//2
        for J in range(0,B+1):
            mf=Mj_formula(A,B,J); mt=Mj((A,B,2),J,m); tested+=1
            if mf is None or mf.denominator!=1 or int(mf)!=mt:
                bad+=1
                if bad<8: print("BAD",(A,B),J,mf,mt)
print(f"interior (b>=2) verify: tested={tested} bad={bad}")

# boundary values j=b+1,b+2 (guarded)
print("--- boundary: M_b, M_{b+1}, M_{b+2}, M_{b+3} for sample shapes (b>=2) ---")
for (A,B) in [(6,4),(8,4),(8,6),(10,6),(7,5),(9,5),(9,7),(12,8)]:
    m=(A+B+2)//2
    row=[]
    for J in range(B,min(B+4,m+1)):
        try: row.append((J,Mj((A,B,2),J,m)))
        except: row.append((J,'ERR'))
    print(f"a={A} b={B} m={m}: {row}")
