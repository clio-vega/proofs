"""Three nested candidates for 'what MS's pairing sees', vs E = Thm 4.1 letters.
   B  = { min L_1^(x) : x in X }          (the actual e~_1 letter)
   U  = union over x in X of L_1^(x)      (every unpaired letter, any x)
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

def L1_of(S,T,x,n):
    return ms_pairing(S,T,x,n)[0]

tot=collections.Counter(); ex=[]
for n in range(3,8):
    c=collections.Counter()
    for S,T,_ in additive_pairs(n):
        X=X_set(S,T,n)
        if not X: continue
        E=clio_letters(S,T,n)
        U=set(); B=set()
        for x in X:
            L=L1_of(S,T,x,n)
            U|=L
            if L:
                pos={r:i for i,r in enumerate(ms_order(x,n))}
                B.add(min(L,key=lambda r:pos[r]))
        c['N']+=1
        c['U_sub_E'] += (U<=E); c['U_eq_E'] += (U==E)
        c['B_sub_U'] += (B<=U)
        if not U<=E and len(ex)<6: ex.append((n,S,T,X,E,U,B))
    print(f"n={n}: N={c['N']}  U subset E: {c['U_sub_E']}  U == E: {c['U_eq_E']}  B subset U: {c['B_sub_U']}")
    for k in c: tot[k]+=c[k]
print("TOTAL",dict(tot))
for n,S,T,X,E,U,B in ex:
    print(f"  U-not-in-E: n={n} S={sorted(S)} T={sorted(T)} X={X} E={sorted(E)} U={sorted(U)} B={sorted(B)}")
