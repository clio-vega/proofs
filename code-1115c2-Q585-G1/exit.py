"""For each changed second journey: (i) is the pre-contact prefix identical (as hexagons) to the solo journey's prefix?
(ii) after the last contact step, is the remaining hexagon sequence identical to a suffix of the solo journey?  (exits cleanly)
(iii) does the journey end while in contact?  Split by BC/BB."""
import pickle
from collections import Counter
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb'))
st = Counter()
for d in D:
    for first, second, alone in (('h1x','h1y','h2y'), ('h2y','h2x','h1x')):
        h2, ha = d[second], d[alone]
        H2 = [s[1] for s in h2]; HA = [s[1] for s in ha]
        if H2 == HA: continue
        solo = set(t for S in HA for t in S)
        contact = [i for i, S in enumerate(H2) if any(t[0] != 'rh' and t not in solo for t in S)]
        f, l = contact[0], contact[-1]
        k = d['kind']
        st[(k, 'prefix identical')] += H2[:f] == HA[:f]
        tail = H2[l+1:]
        st[(k, 'tail is suffix of solo')] += (len(tail) == 0 or HA[len(HA)-len(tail):] == tail)
        st[(k, 'ends in contact')] += (l == len(H2)-1)
        st[(k, 'contact steps contiguous')] += contact == list(range(f, l+1))
        st[(k, 'n')] += 1
for x in sorted(st): print(x, st[x])
