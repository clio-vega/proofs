"""
Q150 structural reduction.

With a = e-1, q = f-e-1, n = a+q+1 = f-1, and sigma, rho not identically zero,
the commutator equations (1B)+(2B) of the Q149 paper are EQUIVALENT to

  (C)   sigma(y[1..a]) rho(y[a+2..n])  =  rho(y[1..q]) sigma(y[q+2..n])   for all y
 (2B*)  sigma(pi 0 pi') sigma(pi' 0 tau) = sigma(pi 1 pi') sigma(pi' 1 tau)
        for 1<=delta<=a, |pi|=|tau|=delta-1, |pi'|=a-delta
        (distinguished letter at position delta on the left factor,
         at position e-delta = a+1-delta on the right factor)

with sigmabar(y) = sigma(y[1..a]) rho(y[a+2..n]).
This module builds the two systems and searches F_p exhaustively.
"""
from itertools import product

def W(k):  return [tuple(w) for w in product((0,1), repeat=k)]

def twobead_pairs(a):
    """list of ((u,v),(u2,v2)) index pairs: sigma[u]sigma[v] = sigma[u2]sigma[v2]"""
    out = []
    e = a+1
    for delta in range(1, a+1):
        for pi in W(delta-1):
            for pip in W(a-delta):
                for tau in W(delta-1):
                    u0 = pi + (0,) + pip ; u1 = pi + (1,) + pip
                    v0 = pip + (0,) + tau; v1 = pip + (1,) + tau
                    assert len(u0)==a and len(v0)==a
                    out.append(((u0,v0),(u1,v1)))
    return out

def compat_pairs(a, q):
    """list of ((uA,vR),(uR,vA)) : sigma[uA]rho[vR] = rho[uR]sigma[vA]"""
    n = a+q+1
    out = []
    for y in W(n):
        lhs = (y[0:a],           y[a+1:n])
        rhs = (y[0:q],           y[q+1:n])
        out.append((lhs, rhs))
    return out

def sigmabar(sig, rho, a, q):
    n = a+q+1
    return {y: sig[y[0:a]]*rho[y[a+1:n]] for y in W(n)}

def solve_Fp(a, q, p):
    """All (sigma,rho) over F_p, both nonzero, satisfying (2B*) and (C)."""
    tb = twobead_pairs(a); cp = compat_pairs(a, q)
    Wa, Wq = W(a), W(q)
    sols = []
    for svals in product(range(p), repeat=len(Wa)):
        sig = dict(zip(Wa, svals))
        if all(v == 0 for v in svals): continue
        if any((sig[x]*sig[y] - sig[z]*sig[w]) % p for (x,y),(z,w) in tb): continue
        # (C) is linear in rho for fixed sigma: brute force rho too (q small)
        for rvals in product(range(p), repeat=len(Wq)):
            if all(v == 0 for v in rvals): continue
            rho = dict(zip(Wq, rvals))
            if any((sig[uA]*rho[vR] - rho[uR]*sig[vA]) % p
                   for (uA,vR),(uR,vA) in cp): continue
            sols.append((sig, rho))
    return sols

# ---------- efficient scan: (C) is LINEAR in rho for fixed sigma ----------
def nullspace_mod(rows, ncols, p):
    M = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(x*inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(M[i][j]-f*M[r][j]) % p for j in range(ncols)]
        piv.append(c); r += 1
        if r == len(M): break
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        v = [0]*ncols; v[fc] = 1
        for i, c in enumerate(piv): v[c] = (-M[i][fc]) % p
        basis.append(v)
    return basis

def scan(a, q, p):
    """All (sigma,rho) over F_p with both nonzero, (2B*) and (C).  Returns
    (n_solutions, n_with_a_zero, list of (support(sigma),support(rho)) pairs)."""
    tb = twobead_pairs(a); cp = compat_pairs(a, q)
    Wa, Wq = W(a), W(q); idx = {v: i for i, v in enumerate(Wq)}
    tot = 0; wz = 0; pats = set()
    for svals in product(range(p), repeat=len(Wa)):
        if not any(svals): continue
        sig = dict(zip(Wa, svals))
        if any((sig[x]*sig[y] - sig[z]*sig[w]) % p for (x,y),(z,w) in tb): continue
        rows = []
        for (uA, vR), (uR, vA) in cp:
            row = [0]*len(Wq)
            row[idx[vR]] = (row[idx[vR]] + sig[uA]) % p
            row[idx[uR]] = (row[idx[uR]] - sig[vA]) % p
            if any(row): rows.append(row)
        for rvals in product(range(p), repeat=len(Wq)):
            if not any(rvals): continue
            if any(sum(r[i]*rvals[i] for i in range(len(Wq))) % p for r in rows): continue
            rho = dict(zip(Wq, rvals))
            tot += 1
            z = (0 in svals) or (0 in rvals)
            if z:
                wz += 1
                pats.add((tuple(u for u in Wa if sig[u]), tuple(v for v in Wq if rho[v])))
    return tot, wz, sorted(pats)

def admissible_supports(a):
    tb = twobead_pairs(a); Wa = W(a); out = []
    for mask in range(1, 1 << len(Wa)):
        S = frozenset(Wa[i] for i in range(len(Wa)) if (mask >> i) & 1)
        if all(((x in S and y in S) == (z in S and w in S)) for (x,y),(z,w) in tb):
            out.append(S)
    return out

def zero_solutions(e, f, p, supports=None):
    """For each PROPER admissible support S, find sigma supported in S satisfying
    (2B*), then solve (C) linearly for rho != 0.  Returns list of (sigma,rho)."""
    a, q = e-1, f-e-1
    tb = twobead_pairs(a); cp = compat_pairs(a, q)
    Wa, Wq = W(a), W(q); idx = {v: i for i, v in enumerate(Wq)}
    if supports is None:
        supports = [S for S in admissible_supports(a) if len(S) < len(Wa)]
    found = []
    for S in supports:
        Sl = sorted(S)
        for vals in product(range(1, p), repeat=len(Sl)):
            sig = {u: 0 for u in Wa}
            sig.update(dict(zip(Sl, vals)))
            if any((sig[x]*sig[y] - sig[z]*sig[w]) % p for (x,y),(z,w) in tb): continue
            rows = []
            for (uA, vR), (uR, vA) in cp:
                row = [0]*len(Wq)
                row[idx[vR]] = (row[idx[vR]] + sig[uA]) % p
                row[idx[uR]] = (row[idx[uR]] - sig[vA]) % p
                if any(row): rows.append(row)
            ns = nullspace_mod(rows, len(Wq), p) if rows else \
                 [[1 if i==j else 0 for i in range(len(Wq))] for j in range(len(Wq))]
            if ns:
                found.append((sig, dict(zip(Wq, ns[0])), S, len(ns)))
    return found

def period_ok(u, e):
    """(2B*) for sigma = c * 1_{u}:  for every period delta of u, u_delta != u_{e-delta}."""
    a = e-1
    for delta in range(1, a+1):
        per = all(u[j] == u[j+delta] for j in range(a-delta))
        if per and u[delta-1] == u[e-delta-1]:
            return False
    return True

def singleton_scan(e, f, p=5):
    """sigma = 1_{u*}; solve (C) linearly for rho != 0.  Returns list of (u*, rho)."""
    a, q = e-1, f-e-1
    cp = compat_pairs(a, q); Wq = W(q); idx = {v: i for i, v in enumerate(Wq)}
    out = []
    for us in W(a):
        if not period_ok(us, e): continue
        sig = {u: (1 if u == us else 0) for u in W(a)}
        rows = []
        for (uA, vR), (uR, vA) in cp:
            row = [0]*len(Wq)
            row[idx[vR]] = (row[idx[vR]] + sig[uA]) % p
            row[idx[uR]] = (row[idx[uR]] - sig[vA]) % p
            if any(row): rows.append(row)
        ns = nullspace_mod(rows, len(Wq), p) if rows else None
        if ns: out.append((us, dict(zip(Wq, ns[0])), len(ns)))
    return out

def tensor(sig_d, d, k):
    """sig_d : {0,1}^{d-1} -> K   |-->   its k-fold block product, rank kd."""
    out = {}
    for y in W(k*d-1):
        v = 1
        for i in range(k):
            v *= sig_d[tuple(y[i*d:(i+1)*d-1])]
        out[y] = v
    return out

def check_system(e, f, sig, rho):
    a, q = e-1, f-e-1
    tb = twobead_pairs(a); cp = compat_pairs(a, q)
    n2 = sum(1 for (x,y),(z,w) in tb if sig[x]*sig[y] != sig[z]*sig[w])
    nC = sum(1 for (uA,vR),(uR,vA) in cp if sig[uA]*rho[vR] != rho[uR]*sig[vA])
    return n2, nC

def twobead_pairs_strong(d):
    """(2B**_d): tau(alpha b alpha') tau(beta b beta') independent of b, with
    the marked letter at position delta0 in the first factor and d-delta0 in the
    second, and alpha,alpha',beta,beta' ALL free and independent.
    This is what (2B*_e) becomes for sigma = tau^{(x)k}, k>=2, at delta = pd+delta0
    with p <= k-2 (the two d-windows are then disjoint)."""
    out = []
    for d0 in range(1, d):
        for al in W(d0-1):
            for alp in W(d-d0-1):
                for be in W(d-d0-1):
                    for bep in W(d0-1):
                        u0 = al+(0,)+alp; u1 = al+(1,)+alp
                        v0 = be+(0,)+bep; v1 = be+(1,)+bep
                        out.append(((u0,v0),(u1,v1)))
    return out

def admissible_supports_gen(d, pairs):
    Wd = W(d-1); out = []
    for mask in range(1, 1 << len(Wd)):
        S = frozenset(Wd[i] for i in range(len(Wd)) if (mask >> i) & 1)
        if all(((x in S and y in S) == (z in S and w in S)) for (x,y),(z,w) in pairs):
            out.append(S)
    return out
