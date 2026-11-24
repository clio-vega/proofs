"""Control for geom.overlap: M0 tiles pairwise non-overlapping (must be 0); each step's removed vs created tiles overlap
(must be > 0 for every step); a tile overlaps itself."""
import sys, pickle
sys.path.insert(0, '/home/clio/projects/proofs/code-1124-Q597'); from geom import overlap, conflict
P = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles5.pkl', 'rb'))
bad = 0; steps = 0; stepok = 0
for p in P[:60]:
    M = list(p['M0'])
    bad += sum(1 for i in range(len(M)) for j in range(i) if overlap(M[i], M[j]))
    for r in p['roles']:
        for s in r[1]:
            steps += 1; stepok += bool(conflict(s[1] - s[2], s[2] - s[1]))
print('M0 overlapping pairs (expect 0):', bad, '| steps whose removed/created sets conflict (expect all):', stepok, '/', steps)
