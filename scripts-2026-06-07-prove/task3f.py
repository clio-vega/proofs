"""
TASK 3f: MO#286705 "generalized binary Krawtchouk" closed form with d=2:
  P_MO(MOm; i, n) = sum_k (-1)^{MOm - k*d} C(n-i, k) C(2i-n, MOm - k*(d+1))   with d=2
Here MOm is the degree-index (we explore = OUR b), and (i,n) are parameters.

We test whether, holding (i,n) fixed, the sequence P_MO(b; i,n) over b matches I_b(OUR-m) ??
But I_b is a polynomial in OUR-m, whereas P_MO is a NUMBER for each (b,i,n).
The instruction: explore whether I_b(OUR-m) as a poly in OUR-m matches their family
in a TRANSPOSED sense -- i.e. is there a dictionary (i,n)<->(b, our-m) so that the
MO number P_MO(b; i,n) reproduces I_b(our-m) when (i,n) are set as affine functions
of our-m?

Concretely: I_b(m0) is an integer for integer m0. We tabulate I_b(m0) for a grid of
(b, m0) and search for affine maps i = a1*m0 + a2*b + a3, n = c1*m0+c2*b+c3 (small int
coeffs) and a global sign/scalar such that P_MO(b; i, n) == I_b(m0) on the grid.
We brute force small integer (i,n) per (b,m0) cell and record which (i,n) reproduce
the value, then look for a consistent linear pattern.
"""
import sympy as sp
from sympy import binomial, Integer
from identify_sweep import I_b, m as M

def P_MO(MOm, i, n, d=2):
    tot = Integer(0)
    for k in range(0, MOm//1 + 2):
        c1 = binomial(n - i, k)
        c2 = binomial(2*i - n, MOm - k*(d+1))
        tot += (-1)**(MOm - k*d) * c1 * c2
    return tot

def main():
    # tabulate I_b(m0)
    grid = {}
    for b in range(2, 11):
        Ib = I_b(b)
        for m0 in range(0, 16):
            grid[(b,m0)] = int(Ib.subs(M, m0))

    # For each (b,m0) find all (i,n) in a box reproducing the value (and sign variants)
    print("Searching (i,n) in box [-20,40] reproducing I_b(m0)=P_MO(b;i,n) ...")
    hits = {}
    for (b,m0), val in grid.items():
        local = []
        for i in range(-20, 41):
            for n in range(-20, 41):
                v = P_MO(b, i, n)
                if v == val:
                    local.append((i,n,+1))
                elif v == -val:
                    local.append((i,n,-1))
        hits[(b,m0)] = local

    # Report cells with small numbers of hits and look for linear (i,n) vs (b,m0)
    # Print a few representative cells
    for (b,m0) in sorted(grid):
        val = grid[(b,m0)]
        h = hits[(b,m0)]
        if val != 0 and 0 < len(h) <= 30:
            print(f" b={b} m0={m0} I={val}: (i,n,sgn) hits = {h[:12]}{' ...' if len(h)>12 else ''}")

    # Try a direct natural dictionary guess: the MO answer there counts a fiber.
    # Common guess: n = 2*m0 (or n related to m0), i related to b. Test n=2*m0, i=m0 etc.
    print("\nTesting specific natural dictionaries P_MO(b; i,n) vs I_b(m0):")
    guesses = [
        ("i=m0, n=2*m0", lambda b,m0:(m0, 2*m0)),
        ("i=m0, n=2*m0+1", lambda b,m0:(m0, 2*m0+1)),
        ("i=m0+b, n=2*m0", lambda b,m0:(m0+b, 2*m0)),
        ("i=b, n=2*m0", lambda b,m0:(b, 2*m0)),
        ("i=m0, n=2*m0-b", lambda b,m0:(m0, 2*m0-b)),
        ("i=2*m0-b, n=2*m0", lambda b,m0:(2*m0-b, 2*m0)),
    ]
    for label, fn in guesses:
        ok=0; tot=0; sgn_consistent=True; first_sgn=None
        for (b,m0),val in grid.items():
            if m0 < b:
                continue
            i,n = fn(b,m0)
            v = P_MO(b,i,n)
            tot+=1
            if v==val: ok+=1
            elif v==-val:
                ok+=1
        print(f"  {label}: matched(±) {ok}/{tot}")

if __name__=="__main__":
    main()
