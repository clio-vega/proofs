"""Verify the PROOF of Theorem M, not just its statement.

The proof says: given y,y' in S and i with y_i > y'_i, put D = {l : y_l < y'_l} and
  Case 1  y_i >= 1, y'_i <= -1 : EVERY l in D works.
  Case 2  y_i >= 1, y'_i >=  0 : every l in D with y'_l >= 1 works, and such an l exists
                                 (else k(y) > k(y'), contradicting k(y)<=j=k(y'));
                                 if k(y) < j every l works.
  Case 3  y_i <= 0 (so y'_i<=-1): every l in D with y_l <= -1 works, and such an l exists
                                 (else k(y) < k(y'), contradicting k(y)=j);
                                 if k(y) < j every l works.
For every (S,y,y',i) encountered: record which case fires, check that the l the PROOF
names is in D, that the exchange it names really lands in S, and that the existence
claim the proof makes is true.  A proof step that is never exercised is reported too."""
import itertools
from collections import Counter
def k(y): return sum(max(-t,0) for t in y)
st=Counter(); bad=[]
for m,lo,hi in [(2,-3,3),(3,-2,2),(4,-2,2)]:
    rngs=[(P,Q) for P in range(lo,hi+1) for Q in range(P,hi+1)]
    for B in itertools.product(rngs,repeat=m):
        for sigma in range(sum(P for P,_ in B), sum(Q for _,Q in B)+1):
            Y=[y for y in itertools.product(*[range(P,Q+1) for (P,Q) in B]) if sum(y)==sigma]
            if len(Y)<2: continue
            for j in sorted({k(y) for y in Y}):
                S=set(y for y in Y if k(y)<=j)
                if len(S)<2: continue
                for y in S:
                    for yp in S:
                        for i in range(m):
                            if y[i]<=yp[i]: continue
                            D=[l for l in range(m) if y[l]<yp[l]]
                            if not D: bad.append(('D empty',B,sigma,j,y,yp,i)); continue
                            # which case?
                            if y[i]>=1 and yp[i]<=-1: case,cand=1,D
                            elif y[i]>=1:              case,cand=2,[l for l in D if yp[l]>=1]
                            else:                      case,cand=3,[l for l in D if y[l]<=-1]
                            st['case%d'%case]+=1
                            if case==2 and k(yp)<j: cand=D; st['case2_slack']+=1   # slack is on the y' side in case 2
                            if case==3 and k(y)<j:  cand=D; st['case3_slack']+=1
                            if not cand:
                                st['case%d_EXISTENCE_FAILED'%case]+=1
                                bad.append(('no candidate',case,B,sigma,j,y,yp,i,D)); continue
                            # the proof claims EVERY named candidate works
                            for l in cand:
                                a=list(y); a[i]-=1; a[l]+=1
                                b=list(yp); b[i]+=1; b[l]-=1
                                ok = tuple(a) in S and tuple(b) in S
                                st['case%d_ok'%case if ok else 'case%d_EXCHANGE_FAILED'%case]+=1
                                if not ok: bad.append(('exchange',case,B,sigma,j,y,yp,i,l))
for key in sorted(st): print(f"  {key} = {st[key]}")
print("  FAILURES:",len(bad))
for b in bad[:5]: print("   ",b)
unexercised=[c for c in (1,2,3) if st['case%d'%c]==0]
print("  unexercised proof cases:", unexercised or "none -- all three cases fire")
