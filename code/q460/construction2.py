"""
Q460 -- the REPAIRED construction and its mechanism check.

Columns 1..4a+2.
  PREFIX rows 1..a  on the even columns  E = {2,4,...,4a}   (2a cols, a arcs)
  SUFFIX rows a+1..2a+1 on  O = {1,3,...,4a+1} u {4a+2}     (2a+2 cols, a+1 arcs)
So N = 2a+1, M = 4a+2, S = (1^{4a+2}).

PROBES (all available in O for every a):
  R_c    = (c-1, c+1)      for every even c in E    [c+1 <= 4a+1, so legal]
  R_inf  = (4a+1, 4a+2)    lies strictly right of every prefix arc

LEMMA (probe table).  For a prefix arc u=(u1,u2) in E and even c in E:
      c <  u1  -> A^T       c = u1 -> C^T
      u1< c <u2-> B         c = u2 -> C        c > u2 -> A
Hence c -> type(u,R_c) determines u1 and u2 (u1 = the unique c giving C^T,
u2 = the unique c giving C), so distinct arcs give distinct probe tables.
And type(u,R_inf) = A for every u.

The order-vector o(type) is injective (divisors.py), so the proportionality
test  "Q -> o(type(P_i,Q_j)) - o(type(P'_i,Q_j)) is constant"  fails whenever
the probe tables differ: compare the separating probe against R_inf, where the
difference is 0.
"""
import itertools, sys

def typ(u, v):
    """tie-free type of (U,V) = (u,v), as the paper's case label."""
    u1, u2 = u; v1, v2 = v
    if u2 < v1: return 'A'
    if v2 < u1: return 'A^T'
    if u1 < v1 and v2 < u2: return 'B'
    if v1 < u1 and u2 < v2: return 'B^T'
    if u1 < v1 < u2 < v2: return 'C'
    if v1 < u1 < v2 < u2: return 'C^T'
    raise ValueError((u, v))

OVEC = {'A':(-1,-1,0), 'A^T':(-1,0,-1), 'B':(-2,0,0), 'B^T':(-2,0,0),
        'C':(-2,-1,1), 'C^T':(-2,1,-1)}

def matchings(cols):
    cols = sorted(cols)
    if not cols:
        yield (); return
    first, rest = cols[0], cols[1:]
    for k in range(len(rest)):
        for m in matchings(rest[:k] + rest[k+1:]):
            yield ((first, rest[k]),) + m

def ordered_matchings(cols):
    return [p for m in matchings(cols) for p in itertools.permutations(m)]

def cols(a):
    E = [2*i for i in range(1, 2*a+1)]                 # 2,4,...,4a
    O = [2*i-1 for i in range(1, 2*a+2)] + [4*a+2]     # 1,3,...,4a+1, 4a+2
    return E, O

print("PART 1 -- probe table lemma: does c -> type(u,R_c) determine u?")
for a in (1,2,3,4,5):
    E, O = cols(a)
    arcs = [(p,r) for p in E for r in E if p < r]
    tables = {}
    ok = True
    for u in arcs:
        tab = []
        for c in E:
            R = (c-1, c+1)
            assert R[0] in O and R[1] in O, (a, c, R)
            tab.append(typ(u, R))
        # check the lemma's predicted shape
        pred = []
        for c in E:
            if   c < u[0]: pred.append('A^T')
            elif c == u[0]: pred.append('C^T')
            elif c < u[1]: pred.append('B')
            elif c == u[1]: pred.append('C')
            else: pred.append('A')
        if tab != pred: ok = False; print("    shape MISMATCH u=%s" % (u,))
        tables[u] = tuple(tab)
    inj = len(set(tables.values())) == len(arcs)
    Rinf = (4*a+1, 4*a+2)
    allA = all(typ(u, Rinf) == 'A' for u in arcs)
    print("  a=%d  |arcs|=%-3d probe-table shape %s  injective %s  R_inf always A %s"
          % (a, len(arcs), "OK" if ok else "FAIL", inj, allA))

print("\nPART 2 -- full separation of all prefix TUPLES by the order-vector test")
for a in (1,2,3):
    E, O = cols(a)
    Ps = ordered_matchings(E)
    suffarcs = [(p,r) for p in O for r in O if p < r]
    tot = good = 0
    for P, Pp in itertools.combinations(Ps, 2):
        tot += 1
        i = next(k for k in range(a) if P[k] != Pp[k])
        u, up = P[i], Pp[i]
        # the proof's prescription: some probe R_c plus R_inf
        d_inf = tuple(s-t for s, t in zip(OVEC[typ(u,(4*a+1,4*a+2))],
                                          OVEC[typ(up,(4*a+1,4*a+2))]))
        sep = False
        for c in E:
            R = (c-1, c+1)
            d = tuple(s-t for s,t in zip(OVEC[typ(u,R)], OVEC[typ(up,R)]))
            if d != d_inf: sep = True; break
        good += sep
    print("  a=%d  %d/%d prefix pairs separated   (|prefixes| = %d = (2a)!/2^a)"
          % (a, good, tot, len(Ps)))

print("\nPART 3 -- NEGATIVE CONTROL: block split (prefix left, suffix right)")
for a in (2,3):
    E = list(range(1, 2*a+1)); O = list(range(2*a+1, 4*a+1))
    Ps = ordered_matchings(E)
    suffarcs = [(p,r) for p in O for r in O if p < r]
    tot = good = 0
    for P, Pp in itertools.combinations(Ps, 2):
        tot += 1
        i = next(k for k in range(a) if P[k] != Pp[k])
        u, up = P[i], Pp[i]
        vals = {tuple(s-t for s,t in zip(OVEC[typ(u,R)], OVEC[typ(up,R)]))
                for R in suffarcs}
        if len(vals) > 1: good += 1
    print("  a=%d  %d/%d separated -> %s" % (a, good, tot,
          "control FIRED (block split separates nothing)" if good==0
          else "control did NOT fire"))
