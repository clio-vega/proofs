import sys; sys.path.insert(0,'/home/clio/projects/proofs/code-1013c2-q442')
from poly import *
from q442 import *
def Ainf_Abd(k,n,w):
    m=2*k; dpU=engineC(k,n,w,False); dpB=engineC(k,n,w,True); Ai={}; Ab={}
    for lam,p in dpU.items():
        if c_minus(lam,m,w)==1: Ai=add(Ai,p)
    for lam,p in dpB.items():
        if c_minus(lam,m,w)==1: Ab=add(Ab,p)
    return Ai,Ab
w=4
for n in range(1,11):
    Ai,Ab=Ainf_Abd(1,n,w)
    d=sub(Ai,Ab)
    al=tuple([1]*n)
    print(f"k=1 w={w} n={n}: A^inf==A {Ai==Ab}; #alpha with A^inf>A = {len(d)}; "
          f"alpha=(1^{n}): A^inf={Ai.get(al,0)} A={Ab.get(al,0)}", flush=True)
print("prediction: first n with A^inf != A is 2(w+1) = %d" % (2*(w+1)))
