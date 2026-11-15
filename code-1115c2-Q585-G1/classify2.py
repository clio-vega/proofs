"""Full crossing word: the sequence of local signatures of ALL steps of the second journey that touch tiles its solo journey never
sees (contact steps), with consecutive repeats compressed.  Also: does the second journey END inside contact (still touching)?"""
import pickle, sys
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from classify import sig_step, cz, D
from collections import Counter, defaultdict
names = {}
def nm(sg):
    if sg not in names: names[sg] = 'K%d' % len(names)
    return names[sg]
out = defaultdict(Counter); exemplars = {}
for idx, d in enumerate(D):
    for first, second, second_alone in (('h1x', 'h1y', 'h2y'), ('h2y', 'h2x', 'h1x')):
        h2, h2a = d[second], d[second_alone]
        if [s[1:] for s in h2] == [s[1:] for s in h2a]: continue
        solo = set(t for s in h2a for t in s[1])
        a = h2a[0]; dv = cz(a[2]) - cz(a[3]); dv /= abs(dv)
        word = []
        for i, s in enumerate(h2):
            new = [t for t in s[1] if t[0] != 'rh' and t not in solo]
            word.append(nm(sig_step(s, set(new), dv)) if new else '.')
        # compress: strip leading/trailing '.' segments, keep pattern
        w = ''.join(x if x == '.' else '[' + x + ']' for x in word).strip('.')
        import re; w = re.sub(r'\.+', '.', w)
        out[w][(d['kind'], d['ok'])] += 1; exemplars.setdefault(w, idx)
print(len(out), 'crossing words')
bc = Counter(); bb = Counter()
for w, c in sorted(out.items(), key=lambda kv: -sum(kv[1].values())):
    print(dict(c), w)
for sg, n in names.items(): print(n, sg)
