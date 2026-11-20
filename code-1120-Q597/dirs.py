import sys, cmath
sys.argv = [sys.argv[0], '5', '1', '1']
exec(open("/home/clio/projects/proofs/code-1119c2-Q597/zonetest.py").read().split("\ndef work(w):")[0])
sys.path.insert(0, '/home/clio/projects/proofs/code-1114c2-Q590'); from common import unit_dir
seen = set()
for idx, (key, b) in enumerate(census_sized()):
    A = b.nest[0]; DV = {'B': A['w'], 'C': A['u']}
    for T in 'BC':
        d = unit_dir(DV[T]); s = (T, DV[T], round(cmath.phase(d)*180/cmath.pi, 3))
        if s not in seen: seen.add(s); print(s)
    if idx > 30: break
print('z0', migrate.G.z0, cmath.phase(migrate.G.z0)*180/cmath.pi, 'z1', cmath.phase(migrate.G.z1)*180/cmath.pi)
