"""Build & pickle primitive integer Q_b for b = 2,3 mod 4, b <= BMAX.
Stored as results/Qb_cache.pkl : dict b -> list of int coeffs (high->low)."""
import sys, time, pickle, os
from _qbcore import Qb_primitive

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    path = "results/Qb_cache.pkl"
    cache = {}
    if os.path.exists(path):
        cache = pickle.load(open(path, "rb"))
    for b in range(2, bmax+1):
        if b % 4 not in (2, 3):
            continue
        if b in cache:
            continue
        t = time.time()
        Q = Qb_primitive(b)
        cache[b] = [int(c) for c in Q.all_coeffs()]
        pickle.dump(cache, open(path, "wb"))   # checkpoint each b
        print(f"b={b:3d} deg={Q.degree():3d} built {time.time()-t:.1f}s", flush=True)
    print(f"cache has {len(cache)} entries -> {path}")

if __name__ == "__main__":
    main()
