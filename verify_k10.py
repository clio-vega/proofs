"""Verify Conjecture 4 for k=10 on S_11 irreps."""
import numpy as np

def partitions(n):
    if n == 0: yield (); return
    def _gen(n, mv):
        if n == 0: yield (); return
        for i in range(min(n, mv), 0, -1):
            for r in _gen(n-i, i): yield (i,) + r
    yield from _gen(n, n)

def standard_tableaux(lam):
    n = sum(lam)
    if n == 0: return [[]]
    results = []
    def fill(t, v):
        if v > n: results.append([list(r) for r in t]); return
        for i in range(len(lam)):
            if len(t[i]) < lam[i]:
                if i == 0 or len(t[i]) < len(t[i-1]):
                    t[i].append(v); fill(t, v+1); t[i].pop()
    fill([[] for _ in range(len(lam))], 1)
    return results

def find_entry(T, val):
    for r, row in enumerate(T):
        for c, v in enumerate(row):
            if v == val: return (r, c)

def axial_distance(T, k):
    r1, c1 = find_entry(T, k)
    r2, c2 = find_entry(T, k+1)
    return (c2 - r2) - (c1 - r1)

def tab_key(T): return tuple(tuple(r) for r in T)

def swap_entries(T, k):
    T_new = [list(r) for r in T]
    r1, c1 = find_entry(T, k); r2, c2 = find_entry(T, k+1)
    T_new[r1][c1] = k+1; T_new[r2][c2] = k
    for r, row in enumerate(T_new):
        for c, v in enumerate(row):
            if c > 0 and v <= row[c-1]: return None
            if r > 0 and c < len(T_new[r-1]) and v <= T_new[r-1][c]: return None
    return T_new

def seminormal_matrix(lam, k):
    tabs = standard_tableaux(lam); d = len(tabs)
    tab_idx = {tab_key(T): i for i, T in enumerate(tabs)}
    M = np.zeros((d, d))
    for T in tabs:
        i = tab_idx[tab_key(T)]; rho = axial_distance(T, k)
        M[i, i] = 1.0 / rho
        T_sw = swap_entries(T, k)
        if T_sw is not None:
            key = tab_key(T_sw)
            if key in tab_idx:
                M[i, tab_idx[key]] = np.sqrt(1.0 - 1.0/rho**2)
    return M

def verify_conjecture4(k):
    n = k + 1
    parts = list(partitions(n))
    print(f"Block k={k}: S_{n} has {len(parts)} partitions")
    
    all_ok = True
    for lam in parts:
        d = len(standard_tableaux(lam))
        if d <= 1: continue
        
        # Staircase through block k-1
        result = np.eye(d)
        for block in range(1, k):
            for j in range(block, 0, -1):
                M = seminormal_matrix(lam, j)
                result = ((np.eye(d) + M) / 2.0) @ result
        
        # Image basis
        U, s, Vt = np.linalg.svd(result, full_matrices=True)
        rank = np.sum(s > 1e-10)
        if rank == 0: continue
        im = U[:, :rank]
        
        # -1(s_k) eigenspace
        Sk = seminormal_matrix(lam, k)
        A = Sk + np.eye(d)
        U2, s2, Vt2 = np.linalg.svd(A)
        null_mask = s2 < 1e-10
        if not np.any(null_mask): continue
        m1 = Vt2[null_mask].T
        
        # Intersection
        combined = np.hstack([im, -m1])
        _, sc, _ = np.linalg.svd(combined)
        inter = max(0, combined.shape[1] - np.sum(sc > 1e-8))
        
        if inter > 0:
            all_ok = False
            print(f"  FAIL: lambda={lam}, dim={d}, rk={rank}, dim_m1={m1.shape[1]}, inter={inter}")
    
    return all_ok

import time

for k in [10]:
    t0 = time.time()
    ok = verify_conjecture4(k)
    t1 = time.time()
    status = "VERIFIED" if ok else "FAILED"
    print(f"  k={k}: {status} ({t1-t0:.1f}s)")
    if ok:
        print(f"  => H-invariant theorem proved for all n <= {k+2}")

