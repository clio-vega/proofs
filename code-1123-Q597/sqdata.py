"""Census corroboration of Lemma sq's mechanism: at the divergence p of every role, the square of A's hexagon H_p is not a
tile of the re-run's hexagon H'_p and overlaps one of its tiles; and it is in W = M0 minus M0.J_o."""
import pickle
exec(open('/home/clio/projects/proofs/code-1122-Q597/geom.py').read())
st = Counter()
for SIZE in (5, 6):
    for P in pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb')):
        M0 = set(P['M0'])
        for (t, A, Bq, Jo) in P['roles']:
            p = next(k for k in range(min(len(A), len(Bq))) if A[k] != Bq[k])
            W = {u for s in Jo for u in s[1] if u in M0}
            H = [it(u) for u in A[p][1]]; H2 = [it(u) for u in Bq[p][1]]
            sq = next(u for u in H if u[0] == 'sq')
            st[('size', SIZE, P['kind'], 'sq not in H\'', sq not in H2, 'overlaps H\'', any(g != sq and overlap(sq, g) for g in H2),
                'sq in W', next(u for u in A[p][1] if u[0] == 'sq') in W)] += 1
for k in sorted(st, key=str): print(st[k], k)
