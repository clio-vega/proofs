"""
Cross-check M_j two ways (build-up Pieri vs SYT-chain removal) for all lambda |- 2m, m <= 8.
The chain model is load-bearing for PROVE, so it MUST agree with the independent build-up.
Also sanity: M_0 = f^lambda, and M_j integer >= 0.
"""
from job_jstar_engine import Mj_buildup, Mj_chain_all, hook_f, partitions

def main():
    bad = 0; checked = 0
    for m in range(1, 9):
        for lam in partitions(2*m):
            Ms_chain = Mj_chain_all(lam)
            for j in range(m + 1):
                a = Mj_buildup(lam, j, m)
                b = Ms_chain[j]
                checked += 1
                if a != b:
                    bad += 1
                    print(f"MISMATCH lam={lam} j={j}: buildup={a} chain={b}")
            # sanity
            if Ms_chain[0] != hook_f(lam):
                bad += 1
                print(f"M_0 != f^lam for lam={lam}: {Ms_chain[0]} vs {hook_f(lam)}")
    print(f"cross-check m<=8: {checked} (lam,j) pairs, {bad} mismatches  ->",
          "PASS" if bad == 0 else "FAIL")

if __name__ == "__main__":
    main()
