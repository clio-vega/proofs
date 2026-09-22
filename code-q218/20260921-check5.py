src=open('/tmp/browse/0921c2-check.py').read()
exec(src[:src.index("for n in range(3,8):")])
from itertools import combinations
def subs_of(n): return [frozenset(c) for r in range(n) for c in combinations(range(n),r)]

def order_list(x,n): return [(x-1-t)%n for t in range(n-1)]   # decreasing order wrt x

def ms_etilde(U,V,x,n):
    """MS Def 3.2(i): u=U left, v=V right.  returns (b, b-t, Unew, Vnew) or None."""
    od=order_list(x,n); pos={z:k for k,z in enumerate(od)}   # smaller pos = larger
    ulist=[z for z in od if z in U]          # decreasing order
    vlist=[z for z in od if z in V]
    usedv=set(); paired=set()
    for b in ulist:
        cands=[a for a in vlist if a not in usedv and pos[a]<pos[b]]   # a > b
        if cands:
            a=cands[-1]        # smallest a>b
            usedv.add(a); paired.add(b)
    L1=[b for b in ulist if b not in paired]
    if not L1: return None
    b=L1[-1]                                  # min of L1 wrt the order
    t=0
    while (b-t-1)%n in U: t+=1
    return b,(b-t)%n

def clio_move(S,T,n):
    out=[]
    for (m,M) in runs(S,n):
        if m in T: continue
        cand=[e for e in cycint(m,M,n) if (e+1)%n not in T]
        if cand: out.append((cand[0],m))
    return out

for n in range(3,8):
    tot=0; agree=0; disagree=[]; msnone=0; clionone=0
    for S in subs_of(n):
        for T in subs_of(n):
            if len(S)<=len(T): continue
            w,L,add=prod(S,T,n)
            if not add: continue
            for x in range(n):
                if x in S or x in T: continue
                tot+=1
                ms=ms_etilde(S,T,x,n)
                cl=clio_move(S,T,n)
                if ms is None: msnone+=1; continue
                if not cl: clionone+=1; continue
                if ms in cl: agree+=1
                else: disagree.append((sorted(S),sorted(T),x,ms,cl))
    print(f"n={n}: (S,T,x) triples {tot} | MS e~ equals SOME Clio (e,m) pair: {agree} | MS gives 0: {msnone} | disagree: {len(disagree)}")
    for d in disagree[:3]: print("    DISAGREE", d)

print("\n--- converse: is every Clio (e,m) realised by MS e~ for some admissible x? ---")
for n in range(3,8):
    tot=0; hit=0; miss=[]
    for S in subs_of(n):
        for T in subs_of(n):
            if len(S)<=len(T): continue
            w,L,add=prod(S,T,n)
            if not add: continue
            xs=[x for x in range(n) if x not in S and x not in T]
            msset={ms_etilde(S,T,x,n) for x in xs}
            for cm in clio_move(S,T,n):
                tot+=1
                if cm in msset: hit+=1
                else: miss.append((sorted(S),sorted(T),cm,sorted(m for m in msset if m)))
    print(f"n={n}: Clio moves {tot}, realised by some MS x: {hit}, not realised: {len(miss)}")
    for m_ in miss[:3]: print("    MISS",m_)
