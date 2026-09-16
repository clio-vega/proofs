"""Enumerate ALL supports Z satisfying the (2B*_d) support condition, by DPLL."""
from twobead import words, relations_2Bstar
import sys

def enumerate_supports(d):
    ws = words(d-1); idx = {w:i for i,w in enumerate(ws)}; n = len(ws)
    rels = []
    seen = set()
    for u0,v0,u1,v1,_ in relations_2Bstar(d):
        key = (idx[u0],idx[v0],idx[u1],idx[v1])
        if key in seen: continue
        seen.add(key); rels.append(key)
    occ = [[] for _ in range(n)]
    for k,(a,b,c,e) in enumerate(rels):
        for v in {a,b,c,e}: occ[v].append(k)
    out = []
    asg = [None]*n

    def propagate(queue):
        """Returns (newly_assigned, ok).  newly is returned even on conflict."""
        newly = []
        while queue:
            v = queue.pop()
            for k in occ[v]:
                a,b,c,e = rels[k]
                A,B,C,D = asg[a],asg[b],asg[c],asg[e]
                # forced assignments
                forced = []
                if A is True and B is True:
                    forced += [(c,True),(e,True)]
                if C is True and D is True:
                    forced += [(a,True),(b,True)]
                if A is False or B is False:
                    if C is True: forced.append((e,False))
                    if D is True: forced.append((c,False))
                if C is False or D is False:
                    if A is True: forced.append((b,False))
                    if B is True: forced.append((a,False))
                for (var,val) in forced:
                    if asg[var] is None:
                        asg[var]=val; newly.append(var); queue.append(var)
                    elif asg[var]!=val:
                        return newly, False
        return newly, True

    def check_all():
        for a,b,c,e in rels:
            if ((asg[a] and asg[b]) != (asg[c] and asg[e])): return False
        return True

    def rec(i):
        while i < n and asg[i] is not None: i += 1
        if i == n:
            if check_all(): out.append(frozenset(ws[j] for j in range(n) if asg[j]))
            return
        for val in (False, True):
            trail = [i]; asg[i]=val
            newly, ok = propagate([i])
            trail += newly
            if ok:
                rec(i+1)
            for v in trail: asg[v]=None
    rec(0)
    return ws, out

if __name__ == '__main__':
    for d in range(3, int(sys.argv[1])+1):
        ws, good = enumerate_supports(d)
        prop = [Z for Z in good if 0 < len(Z) < len(ws)]
        print('d=%d: total admissible %d  (empty 1, full 1, proper nonempty %d)' % (d, len(good), len(prop)))
        bysize = {}
        for Z in prop: bysize[len(Z)] = bysize.get(len(Z),0)+1
        print('      proper by size:', dict(sorted(bysize.items())))
