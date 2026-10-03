"""REFUSAL PANEL for Conjecture N.  Each control's expected outcome is written
   down BEFORE the run, in the `expect` field.  A control that never refuses is
   not a control.

Conjecture N (the statement being tested):
   g >= 0 in Z^m, box B, sigma;  Y = B cap {sum y = sigma};
   w_i(y) = g_i + 1 - (y_i)_+ - (-y_{i+1})_+ >= 1 on Y
   ==>  sum_{y in Y} Trap_{w(y)} is PF_2.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ycoord import centred_sum
from gen import is_pf2
from itertools import product, combinations
from collections import Counter

pos = lambda t: t if t > 0 else 0

def svec(y, pair=1, sign=(1,-1)):
    """s_i = sign-pattern pairing.  pair=1 is the TRUE eq:w pairing."""
    m = len(y)
    a = y[0:]
    return tuple(pos(sign[0]*y[i]) + pos(sign[1]*y[(i+pair) % m]) for i in range(m))

def widths(Y, g, **kw):
    m = len(g)
    ws = [tuple(g[i] + 1 - s[i] for i in range(m)) for s in (svec(y, **kw) for y in Y)]
    return ws if all(min(w) >= 1 for w in ws) else None

def verdict(ws):
    cs = centred_sum(ws)
    if cs is None: return None
    return is_pf2(cs), cs

def run(name, expect, gen_cases, N):
    st = Counter(); wit = []
    for Y, g, kw in gen_cases:
        ws = widths(Y, g, **kw)
        if ws is None: st['zerowidth'] += 1; continue
        v = verdict(ws)
        if v is None: st['parity_gap'] += 1; continue
        ok, cs = v
        st['ok' if ok else 'REFUSED'] += 1
        if not ok and len(wit) < 3: wit.append((g, tuple(Y), sorted(ws), cs))
    tot = st['ok'] + st['REFUSED']
    print(f"  {name:52s} expect={expect:8s} refused {st['REFUSED']}/{tot}   (zerowidth {st['zerowidth']}, parity {st['parity_gap']})")
    for w in wit:
        print(f"      witness g={w[0]} Y={w[1]}")
        print(f"              widths={w[2]} sum={w[3]}")
    return st

# ---------------- case generators ----------------
def boxslices(m, N, R, slack, seed, subset=None, two_planes=False):
    rng = random.Random(seed)
    out = []
    while len(out) < N:
        P = [rng.randint(-R, R) for _ in range(m)]
        Q = [p + rng.randint(0, R) for p in P]
        sig = rng.randint(sum(P), sum(Q))
        rngs = [range(P[i], Q[i]+1) for i in range(m)]
        if two_planes:
            Y = [y for y in product(*rngs) if sum(y) in (sig, sig+1)]
        else:
            Y = [y for y in product(*rngs) if sum(y) == sig]
        if len(Y) < 2 or len(Y) > 400: continue
        if subset == 'random' and len(Y) > 2:
            k = rng.randint(2, len(Y)-1)
            Y = rng.sample(Y, k)
        S = [svec(y) for y in Y]
        maxs = [max(s[i] for s in S) for i in range(m)]
        g = tuple(maxs[i] + rng.randint(0, slack) for i in range(m))
        out.append((Y, g, {}))
    return out

def pairing(m, N, R, slack, seed, pair, sign):
    rng = random.Random(seed); out = []
    while len(out) < N:
        P = [rng.randint(-R, R) for _ in range(m)]
        Q = [p + rng.randint(0, R) for p in P]
        sig = rng.randint(sum(P), sum(Q))
        Y = [y for y in product(*[range(P[i], Q[i]+1) for i in range(m)]) if sum(y) == sig]
        if len(Y) < 2 or len(Y) > 400: continue
        kw = dict(pair=pair, sign=sign)
        S = [svec(y, **kw) for y in Y]
        maxs = [max(s[i] for s in S) for i in range(m)]
        g = tuple(maxs[i] + rng.randint(0, slack) for i in range(m))
        out.append((Y, g, kw))
    return out

if __name__ == '__main__':
    print("=== POSITIVE CONTROL (must NOT refuse): Conjecture N itself ===")
    for m in (3, 4, 5):
        run(f"N, m={m}", "0", boxslices(m, 6000, 5, 3, seed=m), 6000)

    print("\n=== REFUSAL CONTROLS (each must refuse, or it is no control) ===")
    for m in (3, 4):
        run(f"C1 arbitrary SUBSET of the box-slice, m={m}", ">0",
            boxslices(m, 6000, 5, 3, seed=10+m, subset='random'), 6000)
    for m in (3, 4):
        run(f"C2 TWO hyperplanes (sum y in {{s,s+1}}), m={m}", ">0",
            boxslices(m, 6000, 5, 3, seed=20+m, two_planes=True), 6000)
    for m in (4, 5):
        run(f"C3 pairing i,i+2 instead of i,i+1, m={m}", ">0",
            pairing(m, 6000, 5, 3, 30+m, pair=2, sign=(1,-1)), 6000)
    for m in (3, 4):
        run(f"C4 sign (+,+): s_i=(y_i)_+ + (y_{{i+1}})_+, m={m}", ">0",
            pairing(m, 6000, 5, 3, 40+m, pair=1, sign=(1,1)), 6000)
    for m in (3, 4):
        run(f"C5 sign (-,-): s_i=(-y_i)_+ + (-y_{{i+1}})_+, m={m}", ">0",
            pairing(m, 6000, 5, 3, 50+m, pair=1, sign=(-1,-1)), 6000)
