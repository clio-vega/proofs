"""Common prefixes: K=x->B solo, H=x->C solo, G=y->C after K, Gp=y->B after H.  Steps compared as created-tile sets.
Also: y's journeys vs x's journeys translated by -w (one row down)?  usage: python3 prefix.py NSEED"""
import sys
exec(open('/home/clio/projects/proofs/code-1128c2-Q585-SW/swfull.py').read().split("if __name__")[0])
def pre(a, b):
    k = 0
    while k < min(len(a), len(b)) and a[k] == b[k]: k += 1
    return k
def tr(J, v): return [frozenset(tuple(tuple(p+q for p, q in zip(P, v)) for P in vs) for vs in s) for s in J]
st = Counter(); NSEED = int(sys.argv[1])
for k, (key, b0) in enumerate(census()):
    if k >= NSEED: break
    for b, dirs in closure(b0):
        S = set(b.nest[0]['pos']); cs = corners(S); A = b.nest[0]
        for x in cs:
            y = (x[0], x[1]-1)
            if y not in S or y not in corners(S - {x}): continue
            b1 = copy.deepcopy(b); _, K = mv(b1, x, dirs['B']); _, G = mv(b1, y, dirs['C'])
            b2 = copy.deepcopy(b); _, H = mv(b2, x, dirs['C']); _, Gp = mv(b2, y, dirs['B'])
            mw = tuple(-c for c in A['w'])
            st[('pre(K,H)', pre(K, H), 'lenK', len(K), 'lenH', len(H))] += 1
            st[('pre(G,Gp)', pre(G, Gp))] += 1
            st[('pre(G, H shifted -w)', pre(G, tr(H, mw)), 'pre(Gp, K shifted -w)', pre(Gp, tr(K, mw)))] += 1
for k in sorted(st, key=str): print(st[k], k)
