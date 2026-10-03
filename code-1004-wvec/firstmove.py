"""FIRST MOVE (per brief): thm:blind's two certificates. Is the NON-PF2 side M-convex?
If yes, (H2) is dead on arrival."""
import mconv, gen

BAD  = [(1,1),(2,2),(1,5),(2,6)]   # sums to (1,3,4,6,4,3,1), NOT PF2
GOOD = [(1,1),(2,2),(3,3),(4,4)]   # sums to (1,3,6,10,6,3,1), PF2

def conc_sum(fam):
    """sum of concentric convolutions of interval indicators of widths fam[j].
    All centred at 0; returns coefficient list over a common symmetric window."""
    tot = {}
    for w in fam:
        co = gen.conv_intervals(list(w))
        half = (len(co)-1)/2.0
        assert (len(co)-1) % 2 == 0 or True
        for k,c in enumerate(co):
            x = k - half
            tot[x] = tot.get(x,0)+c
    xs = sorted(tot)
    return xs, [tot[x] for x in xs]

print("=== calibrate the sum/PF2 instrument on the two values the node already records ===")
for name, fam, expect_seq, expect_pf2 in [("BAD",BAD,(1,3,4,6,4,3,1),False),
                                          ("GOOD",GOOD,(1,3,6,10,6,3,1),True)]:
    xs, co = conc_sum(fam)
    pf2 = gen.is_pf2(co)
    print(f"  {name}: widths {fam}")
    print(f"     half-widths {[ (sum(w)-len(w))//2 for w in fam]}  (node says 0,1,2,3)")
    print(f"     sum = {tuple(co)}   node says {expect_seq}   match={tuple(co)==expect_seq}")
    print(f"     PF2 = {pf2}   node says {expect_pf2}   match={pf2==expect_pf2}")

print()
print("=== THE TEST: M-convexity of each certificate's width-vector set ===")
for name, fam in [("BAD (non-PF2)",BAD),("GOOD (PF2)",GOOD)]:
    M,wM = mconv.exch_M(fam)
    N,wN = mconv.exch_Mnat(fam)
    print(f"  {name:16s} coord sums {mconv.sums(fam)}")
    print(f"      M-convex   = {M}   witness {wM}")
    print(f"      M^nat-conv = {N}   witness {wN}")
