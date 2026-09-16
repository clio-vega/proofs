from twobead import words, relations_2Bstar, support_ok
from supports import enumerate_supports
from analyse import Cset
import sys

def flip_inv(Z, d, pos):
    """is Z invariant under flipping position pos (1-indexed) of the word?"""
    for w in Z:
        v = list(w); v[pos-1] ^= 1
        if tuple(v) not in Z: return False
    return True

for d in range(3, 8):
    ws, good = enumerate_supports(d)
    prop = [Z for Z in good if 0 < len(Z) < len(ws)]
    # Theorem C'': minimal j with Z <= C_{j,b}; then Z flip-invariant at 1..j-1 and d-j+1..d-1
    ok = True; jcount = {}
    for Z in prop:
        js = [j for j in range(1,d) if 2*j != d and any(Z <= Cset(d,j,b) for b in (0,1))]
        if not js: ok = False; print('   NO C for', Z); continue
        j = min(js); jcount[j] = jcount.get(j,0)+1
        need = list(range(1,j)) + list(range(d-j+1, d))
        for p in need:
            if not flip_inv(Z, d, p):
                ok = False
                print('   FAIL d=%d Z=%s j=%d pos=%d' % (d, sorted(''.join(map(str,x)) for x in Z), j, p))
    print("d=%d: Thm C'' holds on all %d proper supports: %s ; minimal-j distribution %s"
          % (d, len(prop), ok, dict(sorted(jcount.items()))))
    # indicator is always a solution
    rels = list(relations_2Bstar(d))
    bad = 0
    for Z in good:
        f = {w: (1 if w in Z else 0) for w in ws}
        for u0,v0,u1,v1,_ in rels:
            if f[u0]*f[v0] != f[u1]*f[v1]: bad += 1; break
    print('      indicator 1_Z solves (2B*_d) for all %d admissible Z: %s' % (len(good), bad==0))
    sys.stdout.flush()
