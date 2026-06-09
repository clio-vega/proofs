"""
JOB 2c — structure of v2(M_j).  We test:
  (i)  the FULL valuation val(j)=j+2v2(C(m,j)M_j) as a lower-envelope over the
       parity class -- print it for representative ties so the 'box' is visible.
  (ii) a candidate Kummer/Lucas law for v2(M_j).
  (iii) whether the box of J* is already visible in 2v2(M_j)+j (binom removed) or
       needs the binomial carries too.
"""
import math
from job1_tie_census import partitions, M_vector, v2, mn_character

def show(lam,m):
    Ms=M_vector(lam,m)
    print(f"lam={lam} m={m}  M={Ms}")
    par=None
    rows=[]
    for j in range(m+1):
        if Ms[j]==0:
            rows.append((j,'-','-','inf','inf')); continue
        vb=v2(math.comb(m,j)); vM=v2(Ms[j]); val=j+2*vb+2*vM
        rows.append((j,vb,vM,j+2*vM,val))
    # restrict to parity class of the min
    fin=[r for r in rows if r[4]!='inf']
    mu=min(r[4] for r in fin)
    par=[r[0] for r in fin if r[4]==mu][0]%2
    print("  j: vb vM (j+2vM) val   [* = in J*]")
    for (j,vb,vM,jvM,val) in rows:
        if j%2!=par: continue
        star=' *' if val==mu else ''
        print(f"   {j:2d}: {vb}  {vM}   {jvM}    {val}{star}")
    print()

if __name__=="__main__":
    # representative ties of each box type
    for lam,m in [((4,1,1),3),((4,2,2),4),((6,3,1,1,1),6),((8,2,2),6),
                  ((6,6),6),((4,4,2,2),6)]:
        show(lam,m)

    # candidate: is v2(M_j) >= v2(M_0)=v2(f^lam) always? and is the increment structured?
    print("=== v2(M_j) - v2(M_0) increments on a sample (look for additivity) ===")
    import csv
    cnt_mono=0; tot=0
    for m in range(2,11):
        for lam in partitions(2*m):
            Ms=M_vector(lam,m)
            vs=[v2(x) for x in Ms]
            tot+=1
    print(f"(scanned {tot} shapes m<=10)")
