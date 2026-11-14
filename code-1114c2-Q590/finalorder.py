"""Is the FINAL labelled destination tableau independent of the corner-removal order (all linear extensions)?
Labels carried: A-cell (c,r) carries f_+=r+1 when going to B, f_-=alpha'_c - r (depth) when going to C (Lemma flocks).
Also record which orders reproduce R_ii / R_iii (the Purbhoo processing orders)."""
from common import *
from collections import Counter
st = Counter(); witness = []
def runs(b, dvec, want, lab, acc, out):
    I = list(b.nest[0]['pos'])
    if not I:
        out.add(tuple(sorted(acc))); return
    for x in corners(set(I)):
        bx = copy.deepcopy(b); Y, ij, s = mig1(bx, 0, x, dvec); assert Y == want
        runs(bx, dvec, want, lab, acc + [(ij, lab[x])], out)
for key, b0 in census():
    A = b0.nest[0]; uA, wA = A['u'], A['w']; cells = set(A['pos'])
    if len(cells) > 5: st['skipped |alpha|>5'] += 1; continue
    conj = Counter(c for c, r in cells)
    fplus = {(c, r): r+1 for c, r in cells}; fminus = {(c, r): conj[c]-r for c, r in cells}
    for name, dvec, want, lab in (('B', wA, 1, fplus), ('C', uA, 2, fminus)):
        out = set(); runs(copy.deepcopy(b0), dvec, want, lab, [], out)
        st[name+' distinct final labelled tableaux >1'] += len(out) > 1
        st[name+' total outcomes'] += len(out)
        if len(out) > 1 and len(witness) < 2: witness.append((key, name, len(out)))
    st['cells'] += 1
print(dict(st)); print(witness)
