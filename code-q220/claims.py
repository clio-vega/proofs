"""Independent check of the three claims the proof of Theorem A rests on.

For every length-additive (S,T), every legal cut x with L_1 nonempty, write
b = b_x and let [m,M] be the run of S containing b.  Then:
  (C3) b+1 not in T          [local: '(' at b would be closed by ')' at b+1]
  (C2) [m+1,b] subset T      [from minimality of Q(b)]
  (C1) m not in T            [from (C3) + length-additivity via Lemma 4.2 Step 1]
and hence b is the FIRST j in [m,M] with j+1 not in T, i.e. b = e_{[m,M]}.

Also checked: Lemma 4.2 Step 1 itself (m in T  =>  [m,M+1] subset T).
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

cnt=collections.Counter()
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        Tl=set(T)
        # Lemma 4.2 Step 1
        for r in runs(S,n):
            m,M=r[0],r[-1]
            if m in Tl:
                seg=[(m+k)%n for k in range(len(r)+1)]
                cnt['L42S1_tested']+=1
                if not all(z in Tl for z in seg): cnt['L42S1_FAIL']+=1
        for x in X_set(S,T,n):
            res=ms_etilde(S,T,x,n)
            if res is None: continue
            b=res[0]
            r=[rr for rr in runs(S,n) if b in rr][0]
            m,M=r[0],r[-1]
            cnt['fire']+=1
            if (b+1)%n in Tl: cnt['C3_FAIL']+=1
            k=r.index(b)
            if not all(z in Tl for z in r[1:k+1]): cnt['C2_FAIL']+=1
            if m in Tl: cnt['C1_FAIL']+=1
            first=[j for j in r if (j+1)%n not in Tl]
            if not first or first[0]!=b: cnt['EFIRST_FAIL']+=1
            # and the move agrees letter-for-letter with Thm 4.1's
            if (res[1],res[2]) not in clio_moves(S,T,n): cnt['MOVE_FAIL']+=1
print(dict(cnt))

# implementation check on my own side: every Thm 4.1 move really is a factorisation
bad=0; tested=0
for n in range(3,8):
    for S,T,_ in additive_pairs(n):
        w0=window(compose(u_S(S,n),u_S(T,n)),n)
        for Sp,Tp in clio_moves(S,T,n):
            tested+=1
            w1=window(compose(u_S(Sp,n),u_S(Tp,n)),n)
            if w1!=w0 or length(w1,n)!=len(Sp)+len(Tp): bad+=1
print(f"Thm 4.1 moves re-checked as genuine length-additive factorisations: {tested} moves, {bad} failures")
