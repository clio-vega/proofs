src=open('/tmp/browse/0921c2-check.py').read()
exec(src[:src.index("for n in range(3,8):")])
from itertools import combinations
def nilprod(A,B,n):
    w,L,red=word_to_perm(cyc_dec_word(A,n)+cyc_dec_word(B,n),n)
    return w if (red and L==len(A)+len(B)) else None
def iv(m,M,n):
    s=frozenset(cycint(m,M,n)); return None if len(s)>=n else s
for n in range(3,8):
    ok=0;bad=[]
    for i in range(n):
      for j in range(n):
        A=iv(i,j,n)
        if A is None: continue
        for p in range(n):
          for q in range(n):
            B=iv(p,q,n)
            if B is None: continue
            if p in A and j in B:
                if nilprod(B,A,n) is None: ok+=1
                else: bad.append((sorted(B),sorted(A),i,j,p,q))
    print(f"n={n}: Denton ZERO move instances {ok} correct, {len(bad)} FAIL", bad[:2])
