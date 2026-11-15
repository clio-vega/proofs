"""Classify interacting diamonds by the LOCAL configuration at which the second journey first touches the wake of the first.
Wake(J) = tiles (non-rhombus) created by J's rotations that survive to the end of J (Purbhoo §4: tiles displaced by the journey).
Signature of a step = (type, kinds of hexagon tiles, which of them lie in the first journey's wake), with tile centroids
expressed relative to the moving rhombus in the frame where the migration direction points right (rounded floats)."""
import pickle, sys, cmath
sys.path.insert(0,'/home/clio/projects/proofs/code-1112c2-Q585')
from migrate import G
from collections import Counter, defaultdict
D = pickle.load(open('/tmp/c2-diamonds.pkl','rb'))
def cz(vs): return sum(G.phi(P) for P in vs)/len(vs)
def wake(h):
    created = {}; destroyed = set()
    for typ, S, newR, oldR in h:
        for k, vs in S:
            if k != 'rh': destroyed.add((k, vs))
    # tiles after rotation are not logged directly; reconstruct: destroyed tiles minus the rhombus = before-images.
    return destroyed
def sig_step(step, W1, d):
    typ, S, newR, oldR = step
    c0 = cz(oldR); rot = d.conjugate()
    items = []
    for k, vs in S:
        if k == 'rh': continue
        z = (cz(vs) - c0) * rot
        items.append((k, round(z.real, 2), round(z.imag, 2), (k, vs) in W1))
    return (typ, tuple(sorted(items)))
def direction(h):
    typ, S, newR, oldR = h[0]; z = cz(newR) - cz(oldR); return z/abs(z)
out = defaultdict(Counter); first_class = {}
for d in D:
    for first, second, second_alone in (('h1x', 'h1y', 'h2y'), ('h2y', 'h2x', 'h1x')):
        h1, h2, h2a = d[first], d[second], d[second_alone]
        if [s[1:] for s in h2] == [s[1:] for s in h2a]: continue   # this journey not changed by the other
        # post-image wake of the first journey: tiles existing after it that it created
        W1 = set()
        alive = set()
        for typ, S, newR, oldR in h1:
            for k, vs in S: alive.discard((k, vs))
            # we don't log the new tiles; approximate: the tiles destroyed later in h2 that were never in the pre-state
        # simpler: the wake as the set of tiles h2 destroys that h2a (the solo journey) never sees
        solo = set(t for s in h2a for t in s[1])
        dirn = None
        for i, s in enumerate(h2):
            new = [t for t in s[1] if t[0] != 'rh' and t not in solo]
            if new:
                z = cz(s[3]); dd = None
                # direction of migration of the second journey: use the solo journey's net direction
                a = h2a[0]; dv = cz(a[2]) - cz(a[3]); dv /= abs(dv)
                sg = sig_step(s, set(new), dv)
                out[sg][(d['kind'], d['ok'])] += 1
                break
        else:
            out['NO-CONTACT'][(d['kind'], d['ok'])] += 1
print(len(out), 'classes')
for sg, c in sorted(out.items(), key=lambda kv: -sum(kv[1].values())):
    print(dict(c), sg)
