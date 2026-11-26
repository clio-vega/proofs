"""(Mos) at size 6 (REVIEW 11-25 c2 open item): bridge.py logic, unchanged, on EVERY distinct tiling of the
closure states of census chunks lo-hi at size 6 (incl. one-corner states and all one-move images, since the closure is move-closed).
usage: python3 mos6.py NW NCH lo-hi"""
import sys
NW_, NCH_, CH = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sys.argv = ['x', '6', sys.argv[1], sys.argv[2]]
exec(open("/home/clio/projects/proofs/code-1125c2-Q597/bsall.py").read().split("\ndef work(w):")[0])
from fractions import Fraction as Fr
def norm(v): return tuple(int(c) if Fr(c).denominator == 1 else Fr(c) for c in v)
def work(w):
    out = set(); ns = 0; ncorner = Counter()
    for idx, (key, b) in enumerate(census_sized()):
        if idx % NCH != w: continue
        for b, DV in states(b, cap=CAP):
            ns += 1; ncorner[len(corners(set(b.nest[0]['pos'])))] += 1
            out.add(frozenset((k, frozenset(norm(v) for v in vs)) for k, vs in b.tiles.values()))
    return out, ns, ncorner
if __name__ == '__main__':
    T = set(); NS = 0; NC = Counter()
    with Pool(NW_, maxtasksperchild=1) as pool:
        for o, ns, nc in pool.imap_unordered(work, list(range(*map(int, CH.split('-'))))): T |= o; NS += ns; NC.update(nc)
    print('SIZE 6 chunks', CH, 'states', NS, 'corner histogram', dict(NC), 'distinct tilings', len(T), flush=True)
    src = open('/home/clio/projects/proofs/code-1124c2-Q597/bridge.py').read()
    line = "PAIRS = pickle.load(open('/home/clio/projects/proofs/code-1122c2-Q597/roles%d.pkl' % SIZE, 'rb'))"
    assert src.count(line) == 1
    src = src.replace(line, "PAIRS = [{'M0': set(t)} for t in INJECT]")
    sys.argv = ['bridge.py', '6']
    exec(compile(src, 'bridge.py', 'exec'), {'__name__': 'bridge', 'INJECT': sorted(T, key=lambda t: sorted(map(str, t)))})
