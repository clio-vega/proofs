"""X = {} : MS subsection 'Two factor case' (l.1752-1830).

Initial bracketing: bracket i in S with i+1 in T.
Case (3) at index i: a maximal block   u:[b-1..b-t]  v:[b..b-t+1]
   i.e. {b-1,...,b-t} subset S, b not in S,  {b-t+1,...,b} subset T, b-t not in T.
MS: remove the initially bracketed (b-1,b) pair, then proceed with x = b.

Reading A: pair the REDUCED contents S-{b-1}, T-{b} with cut x=b;
           the letter removed is min L_1, and the letter added to T is the bottom
           of its run IN THE ORIGINAL S  (MS's t is defined against content(u)).
Reading B: same, but the bottom of its run in the REDUCED S.
"""
import sys, collections
sys.path.insert(0,'/home/clio/projects/proofs/code-q220')
from affine import *

def case3_bs(S,T,n):
    """all b giving a maximal Case-(3) block."""
    out=[]
    for b in range(n):
        if b in S: continue
        if b not in T: continue
        t=0
        while ((b-t-1)%n in S) and ((b-t)%n in T):
            t+=1
        if t==0: continue
        # maximality is automatic: t is the largest with the block property
        if (b-t)%n in T: continue          # need b-t not in T
        out.append((b,t))
    return out

def ext_move(S,T,b,n,reading='A'):
    So=frozenset(S)-{(b-1)%n}; To=frozenset(T)-{b}
    if b in So or b in To: return None
    L1,R1,_=ms_pairing(So,To,b,n)
    if not L1: return None
    pos={r:i for i,r in enumerate(ms_order(b,n))}
    c=min(L1,key=lambda r:pos[r])
    base = S if reading=='A' else So
    t=0
    while (c-t-1)%n in set(base): t+=1
    m=(c-t)%n
    return c,m

def validate(S,T,Sp,Tp,n):
    if len(Sp)>=n or len(Tp)>=n: return False
    w0=window(compose(u_S(S,n),u_S(T,n)),n)
    w1=window(compose(u_S(Sp,n),u_S(Tp,n)),n)
    return w0==w1 and length(w1,n)==len(Sp)+len(Tp)

tot=collections.Counter(); odd=[]
for n in range(3,8):
    c=collections.Counter()
    for S,T,_ in additive_pairs(n):
        if X_set(S,T,n): continue
        c['pairs']+=1
        E=clio_letters(S,T,n)
        b3=case3_bs(S,T,n)
        c['has_case3'] += bool(b3)
        if not E: c['E_empty']+=1
        for rd in ('A','B'):
            got=set()
            for b,t in b3:
                r=ext_move(S,T,b,n,rd)
                if r: got.add(r)
            letters={x[0] for x in got}
            valid=all(validate(S,T,frozenset(S)-{cc},frozenset(T)|{mm},n) for cc,mm in got)
            c[f'{rd}_valid'] += valid
            c[f'{rd}_sub_E'] += (letters<=E)
            c[f'{rd}_eq_E']  += (letters==E)
            if rd=='A' and not valid and len(odd)<4: odd.append((n,S,T,sorted(got),sorted(E)))
    print(f"n={n} X empty: {c['pairs']} pairs; Case-3 index exists {c['has_case3']}; E empty {c['E_empty']}")
    print(f"     reading A: all moves valid {c['A_valid']}, letters subset E {c['A_sub_E']}, == E {c['A_eq_E']}")
    print(f"     reading B: all moves valid {c['B_valid']}, letters subset E {c['B_sub_E']}, == E {c['B_eq_E']}")
    for k in c: tot[k]+=c[k]
print("TOTAL",dict(tot))
for o in odd: print("  A-invalid: n=%d S=%s T=%s moves=%s E=%s"%(o[0],sorted(o[1]),sorted(o[2]),o[3],o[4]))
