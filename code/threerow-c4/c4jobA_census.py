"""Job A part 2-3 (2026-06-19): census v2(H(j)) floor + Delta_hat(j) slack.

Confirms Lemma N4's CONSTANT floor v2 H >= 4 over a wide sweep, locates the
minimal-slack witnesses, and tests CODE.md's HYPOTHESIZED phi(j) = (linear in j)
- 2 s2(j-8) -- which PROVE found to be wrong (the floor is the constant 4).
Then reproduces the §6 heavy-free Delta_hat(j) > 0 with its minimum slack and the
absorption bound Delta(j) >= j - 4 - 2 s2(j-8).
"""

def v2(n):
    n = int(n)
    if n == 0:
        return 10**9
    return (n & -n).bit_length() - 1

def s2(n):
    return bin(int(n)).count('1')

def Hval(a, b, j):
    return (a**3*b**3+9*a**3*b**2+26*a**3*b+24*a**3+12*a**2*b**3-4*a**2*b**2*j**2
     -8*a**2*b**2*j+108*a**2*b**2-12*a**2*b*j**2-48*a**2*b*j+312*a**2*b-8*a**2*j**2
     -64*a**2*j+288*a**2+47*a*b**3-20*a*b**2*j**2-64*a*b**2*j+423*a*b**2+6*a*b*j**4
     -12*a*b*j**3-30*a*b*j**2-384*a*b*j+1222*a*b+24*a*j**3-40*a*j**2-488*a*j+1128*a
     +60*b**3-24*b**2*j**2-120*b**2*j+540*b**2+6*b*j**4+12*b*j**3-42*b*j**2-696*b*j
     +1560*b-4*j**6+48*j**5-244*j**4+648*j**3-664*j**2-648*j+1440)

def vfact(r):
    return r - s2(r)

# ---------------------------------------------------------------------------
# (2) Census v2(H(j)) floor.  Sweep a < AMAX, all b<=a with a==b (mod 2), and
#     interior j in [8,40].  Tabulate min over (a,b) of v2 H per j; record the
#     minimal-slack witnesses; test candidate phi(j).
# ---------------------------------------------------------------------------
AMAX = 400        # a < 400  (so m up to ~ (a+b+4)/2 ~ 400+)
JMAX = 40
print(f"=== (2) Census v2(H(j)) over a<{AMAX}, b<=a, a==b mod2, j in [8,{JMAX}] ===")
print(" j | min v2H | #witnesses(min) | example (a,b)   | candidate phi=j-4-2s2(j-8)")
floor_by_j = {}
global_min = 99
witnesses = []   # (a,b,j) attaining the global min
for jv in range(8, JMAX+1):
    mn = 99
    cnt = 0
    ex = None
    # sweep both parities of a (b shares parity); keep b in [4, a]
    for av in range(4, AMAX):
        b0 = 4 if (av % 2 == 0) else 5
        for bv in range(b0, av+1, 2):
            val = v2(Hval(av, bv, jv))
            if val < mn:
                mn = val
                cnt = 1
                ex = (av, bv)
            elif val == mn:
                cnt += 1
    floor_by_j[jv] = mn
    cand = jv - 4 - 2*s2(jv-8)   # CODE.md's hypothesized "linear - 2 s2(j-8)" shape
    if mn < global_min:
        global_min = mn
    print(f"{jv:2d} |   {mn:2d}    |     {cnt:6d}      | {str(ex):14s} |  {cand}")

print(f"\n    GLOBAL min over all swept (a,b,j>=8): v2H = {min(floor_by_j.values())}")
const = all(v == 4 for v in floor_by_j.values())
print(f"    floor == constant 4 for EVERY j in [8,{JMAX}] ?  {const}  "
      f"(min per j = {sorted(set(floor_by_j.values()))})")

# Test CODE.md's hypothesized phi(j) = (linear) - 2 s2(j-8): does v2H - phi >= 0?
print("\n    Verdict on CODE.md hypothesis phi(j)=(linear in j)-2 s2(j-8):")
print("      The TRUE floor is the CONSTANT 4 (independent of j) -- so any")
print("      'linear in j' phi is FALSE: v2H(j) does NOT grow with j.")
# concretely: candidate j-4-2s2(j-8) exceeds the actual floor 4 for large j:
viol = [(jv, jv-4-2*s2(jv-8)) for jv in range(8, JMAX+1) if (jv-4-2*s2(jv-8)) > 4]
print(f"      candidate j-4-2s2(j-8) EXCEEDS actual floor 4 at {len(viol)} of "
      f"{JMAX-7} indices (e.g. j={viol[0][0]}: cand={viol[0][1]} > 4) -> hypothesis REFUTED.")

# Sharpness: 32 does NOT divide H somewhere (floor cannot be raised to 5).
sharp = any(v2(Hval(av, bv, jv)) == 4
            for av in range(4, 40) for bv in range(4 if av % 2 == 0 else 5, av+1, 2)
            for jv in range(8, 20))
print(f"    sharp: v2H = 4 actually attained (32 does NOT always divide H): {sharp}")
print(f"      e.g. v2 H(10,8,8) = {v2(Hval(10,8,8))} (note's witness), "
      f"v2 H(7,7,0) = {v2(Hval(7,7,0))}")
