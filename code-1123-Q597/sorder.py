"""Removal orders of J1, J2 on S = W1 cap W2 (M0 tiles removed by both). usage: python3 sorder.py SIZE"""
import sys, pickle
from collections import Counter
SIZE=int(sys.argv[1])
PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))
st=Counter(); ex=0
def remtimes(J,M0):
    r={}
    for k,(lab,rem,cre) in enumerate(J):
        for u in rem:
            if u in M0: assert u not in r; r[u]=k
    return r
for P in PAIRS:
    M0=set(P['M0'])
    J1=P['roles'][0][1]; J2=P['roles'][1][1]   # A of role x is J_x; A of role y is J_y
    r1=remtimes(J1,M0); r2=remtimes(J2,M0)
    S=set(r1)&set(r2)
    p1=min(r1[u] for u in S); p2=min(r2[u] for u in S)
    # orders as sequences of time-pairs, sorted by J1 time
    seq=sorted((r1[u],r2[u],u[0]) for u in S)
    rel=tuple((a-p1,b-p2,k) for a,b,k in seq)
    # control: c0 - p
    c01=min(r2[u] for u in S if r1[u]==p1); c02=min(r1[u] for u in S if r2[u]==p2)
    st[('ctl c0(1)-p2',c01-p2,'c0(2)-p1',c02-p1)]+=1
    st[(P['kind'],'|S|',len(S))]+=1
    st[(P['kind'],'rel',rel)]+=1
for k in sorted(st,key=str): print(st[k],k)
