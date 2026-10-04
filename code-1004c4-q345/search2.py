import reconstruct as R
from itertools import permutations
from schub import schubert_table
IK = ["1..a","1..b-1","a..b-1","a..a","b-1..b-1"]
CD = ["free","weak","strict","weak-dec","strict-dec","bjs","bjs-rev"]
KK = ["a","b","len"]
cands=[]
for ikind in IK:
    for cond in CD:
        for kkind in (KK if cond.startswith("bjs") else ["none"]):
            for adj in (False,True):
                cands.append((ikind,cond,kkind,adj))
print("rules in search space:", len(cands), flush=True)
for n in (3,4):
    tab = schubert_table(n,n)
    targets = {v: R.schub_dict(tab[v],n) for v in permutations(range(1,n+1))}
    keep=[]
    for c in cands:
        if all(R.gen(v,n,c[0],c[1],c[2],c[3])==targets[v] for v in targets):
            keep.append(c)
    print(f"n={n}: {len(keep)} rules reproduce S_v for every permutation:", keep, flush=True)
    cands=keep
    if not cands: break
