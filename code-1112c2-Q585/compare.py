"""Q585: Purbhoo's Cor 3.4 swap M vs the Henriques-Kamnitzer crystal commutor C, on all c>=2 cells
(gl_4, |lam|,|mu| <= 5; the 67 cells of PROVE c1).  Two independent tableau <-> highest-weight bijections.

beta_RSK(b; lam, mu): b⊗u_mu hw in B(lam)⊗B(mu) (word read(b)+read(u_mu)).  Row-insert read(b) (gives b),
   then the letters of read(u_mu); record each new box with the inserted LETTER -> LR tableau on nu/lam, content mu.
beta_CM(b; lam, mu): content-matrix transpose: row r of nu/mu holds m_{k,r} copies of k, where m_{k,r} =
   #(r+1)'s in row k of b -> LR tableau on nu/mu, content lam.
"""
import sys, time
sys.path.insert(0, '/home/clio/projects/proofs/code-1112-Q585')
from commutor import LR, C1, C2, parts, u, rows_read
from crystal import rsk
from flock import cor34, is_lr, content
NC = 4

def french(rows_by_r):   # dict (r,c)->v  ->  dict (col,row)->v
    return {(c, r): v for (r, c), v in rows_by_r.items()}
def beta_rsk(b, lam, mu):
    w = rows_read([list(r) for r in b]); P, Q = rsk(w); assert [len(r) for r in P] == list(lam)
    um = rows_read(u(mu)); P = [list(r) for r in P]; T = {}
    for x in um:
        r = 0; lab = x
        while True:
            if r == len(P): P.append([x]); T[(r, 0)] = lab; break
            row = P[r]; j = next((j for j, y in enumerate(row) if y > x), None)
            if j is None: row.append(x); T[(r, len(row)-1)] = lab; break
            row[j], x0 = x, row[j]; x = x0; r += 1
    return french(T)
def beta_cm(b, lam, mu, nu):
    T = {}
    for r in range(len(nu)):
        vals = sorted(k+1 for k, row in enumerate(b) for y in row if y == r+1)
        m = mu[r] if r < len(mu) else 0
        assert len(vals) == nu[r]-m
        for i, v in enumerate(vals): T[(m+i, r)] = v
    return T
def box(shape_outer, shape_inner):
    return {(c, r) for r, l in enumerate(shape_outer) for c in range(shape_inner[r] if r < len(shape_inner) else 0, l)}
def trim(p): return tuple(x for x in p if x)

cache = {}
def M(T, outer):
    key = (tuple(sorted(T.items())), outer)
    if key not in cache:
        k0 = max(outer[0], 1)+1; k1 = len(outer)+1
        g, info, _ = cor34(T, outer, k0, k1)
        cache[key] = g
    return cache[key]

if __name__ == '__main__':
    stats = dict(cells=0, elems=0, rsk_agree=0, cm_agree=0, lr_bad=0, shape_bad=0, invol_bad=0,
                 rsk_cells_agree=0, cm_cells_agree=0, lex_differs_M=0)
    t0 = time.time(); disagree = []
    for a in range(1, 6):
        for bs_ in range(1, 6):
            for lam in parts(a, NC):
                for mu in parts(bs_, NC):
                    L1 = LR(lam, mu, NC); L2 = LR(mu, lam, NC)
                    for nu, bs in L1.items():
                        if len(bs) < 2: continue
                        nu_t = trim(nu); stats['cells'] += 1
                        ra = ca = True; Mimg = {}
                        for b in bs:
                            stats['elems'] += 1
                            cb = C1(b, lam, mu, NC)[0]
                            assert cb == C2(b, lam, mu, NC)[0] and cb in L2[nu]
                            # bijection 2
                            T2 = beta_cm(b, lam, mu, nu_t); assert is_lr(T2) and content(T2) == trim(lam)
                            out2 = M(T2, nu_t)
                            want2 = beta_cm(cb, mu, lam, nu_t); Mimg[b] = out2
                            if not is_lr(out2): stats['lr_bad'] += 1
                            if set(out2) != box(nu_t, lam): stats['shape_bad'] += 1
                            if M(out2, nu_t) != T2: stats['invol_bad'] += 1
                            ok2 = (out2 == want2); stats['cm_agree'] += ok2; ca &= ok2
                            if not ok2: disagree.append(('cm', lam, mu, nu_t, b, cb))
                        stats['rsk_cells_agree'] += ra; stats['cm_cells_agree'] += ca
                        src = sorted(bs); tgt = sorted(L2[nu])
                        lexmap = {b: beta_cm(tgt[src.index(b)], mu, lam, nu_t) for b in bs}
                        stats['lex_differs_M'] += any(lexmap[b] != Mimg[b] for b in bs)
                        stats['lex_differs_C'] = stats.get('lex_differs_C', 0) + any(tgt[src.index(b)] != C1(b, lam, mu, NC)[0] for b in bs)
    print(stats, 'time %.0fs' % (time.time()-t0))
    for d in disagree[:6]: print('DISAGREE', d)
