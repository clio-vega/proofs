"""PROVE 2026-11-27 c2, Part 2 Step 0: the dictionary 'J_x = Purbhoo migration of the accessible singleton flock {x}'.
On the uncapped size-5 closure (all states, all nest-A corners, both targets) check per journey:
 P1 start cell is a removable corner of nest A (accessible singleton flock);
 P2 the journey lands in the intended nest (w -> B = nest 1, u -> C = nest 2);
 P3 every step moves the rhombus centroid weakly right in the migration frame (x = Re(phi(P) * conj(d)));
 P4 ties: steps where the engine saw >1 candidate hexagon (Purbhoo's 'unique minimal hexagon' fails there);
 P5 every reached board is a Purbhoo mosaic: oriented-edge manifold check + area = area(H) + every nest an ideal;
 P6 lesion: pop a NON-corner nest-A cell, journey toward B; is the result a mosaic? does the psi read-back rebuild it?
usage: python3 dict_falsifier.py NW NCH"""
import sys, copy
from collections import Counter
from fractions import Fraction as Fr
from multiprocessing import Pool
sys.path.insert(0, '/home/clio/projects/reviews/code-1125c2-BS')
sys.path.insert(0, '/home/clio/projects/proofs/code-1112c2-Q585')
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590')
import migrate, mosaic as MZ
from indep_bs_census import census, closure, state_key
from common import corners, tiles_of
from flock import unit_dir
G = migrate.G
U = {(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)}
def psi(P): return (P[0] + P[2], P[1] + P[3])

_step = migrate.Board.step

import inspect, textwrap
LES = [0]
import indep_bs_census as IBC
_s = textwrap.dedent(inspect.getsource(IBC._engine_step))
assert _s.count("k = ks[0]; newR = None") == 1
_s = _s.replace("k = ks[0]; newR = None", "k = ks[0]; newR = None\n    if LES[0] and typ == 'b': k = 6 - k; LES[0] = 0")
_ns = dict(vars(migrate)); _ns['LES'] = LES; exec(_s, _ns); _lesion_step = _ns['step']
LOG = []
def logged(self, R, d):
    t0 = migrate.Global_ties[0]; rot = d.conjugate()
    xb = sum((G.phi(P) * rot).real for P in self.tiles[R][1]) / 4
    R2 = (_lesion_step if LES[0] else _step)(self, R, d)
    if R2 is not None:
        xa = sum((G.phi(P) * rot).real for P in self.tiles[R2][1]) / 4
        LOG.append((self.last, xa - xb, migrate.Global_ties[0] - t0))
    return R2
migrate.Board.step = logged

def area(vs):
    z = [G.phi(P) for P in vs]; s = 0
    for i in range(len(z)): s += (z[i].conjugate() * z[(i + 1) % len(z)]).imag
    return s / 2

def is_mosaic(b, Harea):
    """oriented half-edges: every edge used by <=2 tiles, twice only with opposite orientation; all tiles ccw;
    total area = area(H); boundary edges form one cycle; every nest an ideal."""
    he = Counter(); tot = 0
    for kind, vs in b.tiles.values():
        a = area(vs)
        if a <= 1e-9: vs = vs[::-1]; a = -a
        tot += a
        for i in range(len(vs)): he[(vs[i], vs[(i + 1) % len(vs)])] += 1
    if any(c > 1 for c in he.values()): return False, 'half-edge reused'
    bd = [e for e in he if (e[1], e[0]) not in he]
    nxt = {}
    for p, q in bd:
        if p in nxt: return False, 'boundary not a simple cycle'
        nxt[p] = q
    p0 = bd[0][0]; p = p0; k = 0
    while True:
        p = nxt[p]; k += 1
        if p == p0: break
        if k > len(bd): return False, 'boundary walk'
    if k != len(bd): return False, 'boundary has >1 component'
    if abs(tot - Harea) > 1e-6: return False, 'area'
    if not all(migrate.Board.is_ideal(set(info['pos'])) for info in b.nest): return False, 'nest not ideal'
    inn = [t for info in b.nest for t in info['pos'].values()]
    if sorted(inn) != sorted(t for t, (k, _) in b.tiles.items() if k == 'rh'): return False, 'rhombus outside the nests'
    return True, 'ok'

def bridge_ok(tileset):
    """psi read-back (logic of code-1124c2-Q597/bridge.py, copied): labels off Z^4 differences, KT pieces, rebuild."""
    M0 = tileset
    allp = {psi(Q) for t in M0 if t[0] != 'rh' for Q in t[1]}
    a0 = min(p[0] for p in allp); b0 = min(p[1] for p in allp); n = max(p[0] + p[1] for p in allp) - a0 - b0
    sh = lambda p: (p[0] - a0, p[1] - b0)
    lab = {}; okl = True
    for t in M0:
        if t[0] == 'rh': continue
        V = list(t[1])
        for x in V:
            for y in V:
                D = tuple(bb - aa for aa, bb in zip(x, y))
                if (D[2:] == (0, 0) and D[:2] in U): l = 0
                elif (D[:2] == (0, 0) and D[2:] in U): l = 1
                else: continue
                e = MZ.E(sh(psi(x)), sh(psi(y)))
                if lab.get(e, l) != l: okl = False
                lab[e] = l
        if t[0] == 'sq':
            Vs = [sh(psi(x)) for x in V]
            for x in Vs:
                for y in Vs:
                    if x < y and (y[0] - x[0], y[1] - x[1]) in U and MZ.E(x, y) not in lab: lab[MZ.E(x, y)] = 2
    if not okl: return False
    if not all(tuple(lab.get(MZ.E(c[i], c[(i + 1) % 3])) for i in range(3)) in MZ.OK for c in MZ.tri_cells(n)): return False
    try: M = MZ.build_mosaic(n, lab)
    except Exception: return False
    built = {(k, frozenset(v)) for k, v in M['tiles']} | {(k, frozenset(v)) for nest in M['nests'] for k, v in nest}
    Pm = [Q for t in M0 for Q in t[1] if sh(psi(Q)) == (0, 0) and t[0] != 'rh'][0]
    T0 = M['L'][(0, 0)]; d = tuple(a - b for a, b in zip(Pm, T0))
    built = {(k, frozenset(tuple(a + b for a, b in zip(Q, d)) for Q in v)) for k, v in built}
    return built == set(M0)

def tileset(b): return {(k, frozenset(vs)) for k, vs in b.tiles.values()}

def work(args):
    w, NCH = args; st = Counter()
    for i, (key, b0) in enumerate(census()):
        if i % NCH != w: continue
        Harea = sum(area(vs) if area(vs) > 0 else -area(vs) for _, vs in b0.tiles.values())
        for b, dirs in closure(b0):
            st['states'] += 1
            ok, why = is_mosaic(b, Harea); st[('P5 state is a Purbhoo mosaic', why)] += 1
            S = set(b.nest[0]['pos']); cs = corners(S)
            for x in cs:
                st[('P1 start is removable corner', b.nest[0]['pos'] and migrate.Board.is_ideal(S - {x}))] += 1
                for T, X in (('B', 1), ('C', 2)):
                    nb = copy.deepcopy(b); LOG.clear()
                    t = nb.nest[0]['pos'].pop(x)
                    R, Y, ij, steps = nb.journey(t, 0, unit_dir(dirs[T]))
                    st['journeys'] += 1; st[('P2 lands in intended nest', Y == X)] += 1
                    for kind, dx, ties in LOG:
                        st['steps'] += 1; st[('P3 centroid weakly right', dx > -1e-9)] += 1
                        if abs(dx) < 1e-9: st['P3 strictly vertical step'] += 1
                        if ties: st['P4 tie steps'] += 1
                    if any(t_ for _, _, t_ in LOG): st['P4 journeys with a tie'] += 1
                    ok2, why2 = is_mosaic(nb, Harea); st[('P5 result is a Purbhoo mosaic', why2)] += 1
            # lesion: non-corner cells, toward B (one per state is enough to see it, take all)
            for x in sorted(S - set(cs)):
                nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x)
                try:
                    R, Y, ij, steps = nb.journey(t, 0, unit_dir(dirs['B']))
                except Exception as ex:
                    st[('P6 lesion journey', 'raised ' + str(ex)[:30])] += 1; continue
                st[('P6 lesion journey', 'terminated')] += 1
                ok3, why3 = is_mosaic(nb, Harea); st[('P6 lesion result mosaic?', why3)] += 1
                st[('P6 lesion psi-readback rebuilds board', bridge_ok(tileset(nb)))] += 1
            # lesion L2: first type-b step of the journey uses the OTHER 120-degree rotation (violates the edge rule)
            for x in cs[:1]:
                nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); LES[0] = 1
                try:
                    R, Y, ij, steps = nb.journey(t, 0, unit_dir(dirs['B']))
                except Exception as ex:
                    st[('L2 journey', 'raised ' + str(ex)[:30])] += 1; LES[0] = 0; continue
                if LES[0]: st[('L2 journey', 'no type-b step, lesion did not fire')] += 1; LES[0] = 0; continue
                st[('L2 journey', 'terminated, Y=%d' % Y)] += 1
                ok4, why4 = is_mosaic(nb, Harea); st[('L2 result mosaic?', why4)] += 1
                st[('L2 psi-readback rebuilds board', bridge_ok(tileset(nb)))] += 1
            # lesion L3: stop the journey one step early (the rhombus is never packed)
            for x in cs[:1]:
                nb = copy.deepcopy(b); t = nb.nest[0]['pos'].pop(x); LOG.clear()
                ref = copy.deepcopy(nb); ref.journey(ref.nest[0]['pos'].get(x, t) if False else t, 0, unit_dir(dirs['B']))
                s = len(LOG)
                R = t
                for _ in range(s - 1): R = nb.step(R, unit_dir(dirs['B']))
                ok5, why5 = is_mosaic(nb, Harea); st[('L3 early stop: result mosaic?', why5)] += 1
                st[('L3 early stop: psi-readback rebuilds board', bridge_ok(tileset(nb)))] += 1
            if st['states'] % 50 == 1:   # control for the readback instrument on unlesioned states
                st[('P6 control psi-readback on state', bridge_ok(tileset(b)))] += 1
    return st

if __name__ == '__main__':
    NW, NCH = int(sys.argv[1]), int(sys.argv[2]); tot = Counter()
    with Pool(NW, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(work, [(w, NCH) for w in range(NCH)]): tot.update(o)
    for k in sorted(tot, key=str): print(tot[k], k)
