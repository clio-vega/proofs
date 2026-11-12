"""Bundle crossing B_N: SE-bundle (word s) x SW-bundle (word t), every SE line crosses every SW line once
with the validated Schur-point diamond D0 (reading R2 / 'top,left=SW,right=SE').
D0: (UL=SE_in, UR=SW_in) -> {(LL=SW_out, LR=SE_out)}."""
from itertools import product
from collections import Counter
DELTA={(0,0,0),(1,0,2),(2,1,0),(1,1,1),(0,2,1)}
NAB={(0,0,0),(0,2,1),(2,1,0),(1,0,2),(1,1,1)}
D0={}
for ul,ur,ll,lr in product((0,1,2),repeat=4):
    c=sum((x,ul,ur) in DELTA and (x,lr,ll) in NAB for x in (0,1,2))
    if c: D0.setdefault((ul,ur),[]).append(((ll,lr),c))
def B(s,t):
    """returns Counter of (SE_out word, SW_out word) with multiplicity"""
    # process SE line i=0..N-1 sequentially; it crosses SW lines j=0..M-1 in order
    states=Counter({(tuple(),tuple(t)):1})
    for a in s:
        new=Counter()
        for (seout,sw),mult in states.items():
            # thread SE label a through all SW lines
            cur=Counter({(a,()):mult})
            for b in sw:
                nxt=Counter()
                for (x,swn),mm in cur.items():
                    for (ll,lr),c in D0.get((x,b),[]):
                        nxt[(lr,swn+(ll,))]+=mm*c
                cur=nxt
            for (x,swn),mm in cur.items(): new[(seout+(x,),swn)]+=mm
        states=new
    return states
if __name__=='__main__':
    print('D0 table:',D0)
    for N in (2,3,4):
        tot=killed=single=swapped=0
        for s in product((0,1),repeat=N):
            for t in product((0,1),repeat=N):
                if sum(s)!=sum(t): continue   # same Grassmannian Gr(k,N) for both
                out=B(s,t); tot+=1
                if not out: killed+=1; continue
                if len(out)==1 and sum(out.values())==1:
                    single+=1
                    (so,to),=out
                    swapped+= (so,to)==(t,s) or (so,to)==(t[::-1],s[::-1])
        print(f'N={N}: pairs {tot}, B kills {killed}, single-output {single}, of which = swap {swapped}')
