"""PROVE 2026-11-28, Lemma Cand (finite): every cluster accepted by the engine's hexagon test (hexenum.py's exhaustive
list, 144 clusters) is R + the FULL fan at one obtuse corner of R: the three other tiles all contain that corner v
and their angles at v sum to 360 - 150 = 210.  Reuses hexenum.py's enumeration verbatim (exec), adds the check."""
import cmath, math
from collections import Counter
src = open('/home/clio/projects/proofs/code-1127c2-Q597/hexenum.py').read()
src = src.split("print('accepted clusters")[0]
src = src.replace("    acc += 1\n", "    acc += 1; ACC.append(cl)\n")
ACC = []
exec(src)
def angle_at(cyc, v):
    i = cyc.index(v); a = ph(cyc[i - 1]) - ph(v); b = ph(cyc[(i + 1) % len(cyc)]) - ph(v)
    return abs(math.degrees(cmath.phase(b / a)))
st = Counter()
for cl in ACC:
    R = [c for c in cl if kind(c) == 'rh'][0]; others = [c for c in cl if c is not R]
    obt = [v for v in R if abs(angle_at(R, v) - 150) < 1e-6]
    hits = [v for v in obt if all(v in c for c in others) and abs(sum(angle_at(c, v) for c in others) - 210) < 1e-6]
    st[('accepted cluster = R + full fan at an obtuse corner; #such corners', len(hits))] += 1
print('accepted clusters:', len(ACC))
for k, v in sorted(st.items(), key=str): print(v, k)
