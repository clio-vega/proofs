# c=5 interior PROVED — `(a,b,5)` is now a COMPLETE theorem

**2026-06-19, prove session.** Hi Robin.

The PROVE target was the interior even-|J*| theorem for three-row `(a,b,5)`. **Done, with no
residual.** With the c=5 boundary (closed 06-17), the three-row family `(a,b,5)` is now genuinely
complete — interior *and* boundary — the **first odd-c family complete above c=3**, joining c=1,2,3,4.

Proof: `projects/proofs/2026-06-19-c5-interior-Jstar-even.md`. Code + verification:
`projects/code/threerow-c5/`.

## The headline is a surprise

PROVE.md (and I) expected c=5 to be a two-generator monster like c=3, which has the full box
`|J*|=4`. **It isn't.** c=5, though odd, behaves like c=4: **only the generator 2 fires, `|J*| ≤ 2`
always.** The offset is set by the parity of `a`:

- `a` even (`b` odd): `j₀ = 0`, `J* = {0,2}` iff `(a,b) ≡ (0,1) mod 4`, else `{0}`.
- `a` odd (`b` even): `j₀ = 3` (forced descent Δ(1),Δ(2),Δ(3) = −1,−2,−3, exact), `J* = {3,5}` iff
  `(a,b) ≡ (3,0) mod 4`, else `{3}`.

End-to-end Murnaghan–Nakayama check: 561 box-interior shapes m≤45, 128 ties, **0 mismatch**.

The reason (§9 of the writeup): for c≥4 the tip `(2c)!C(j,2c)` first turns nonzero at `j=2c`, which by
then is already in the "deep régime" where absorption + the heavy floor force Δ>0. So the would-be
generators 4,6,8,10 never get a low-content window to fire. c=3's `|J*|=4` is special to small c (its
tip nonzero point 2c=6 sits *below* the deep threshold).

## For your & Rick's FREE/RIGID note — the new data point

The forward question was: does the heavy quotient `H_c` carry a *constant* 2-adic floor `2^{β'(c)}`
(no a,b-content)? **Yes.** New witness: **β'(5) = 3** (`8 ∣ H_5`, sharp at (a,b,j)=(3,0,2)).

Combined with β'(4) = 4 (`16 ∣ H` from the c=4 Number Lemma), this tells us:
- `β'(c)` is **not monotonic** (4 → 3), and
- it is **unrelated** to the rigid floor `β(c) = (c−1)+v₂((c−1)!) = 2(c−1)−s₂(c−1)` of the single-
  binomial Number Lemma NL_c (which gives β(5)=6). The heavy *product* has its own, smaller, floor.

My reading of the 4→3 drop: it tracks the parity régime. c=4 forces `a ≡ b (mod 2)`, so both even-side
factors of `H(0)` are available; c=5 forces `a,b` opposite, costing one guaranteed factor of 2. If
that's the right lens, I'd predict `β'(c)` alternates with the parity constraint rather than growing —
worth testing at c=6,7.

## What I'd do next (your call)

1. **LEAN c=5** — `8 ∣ H_5` is decide-friendly (residue check mod 16), exactly like the c=4 `16∣H`
   that's already formalised. The whole c=5 interior reduces to a handful of finite residue checks.
2. **General-c single-generator theorem** — the diagnosis suggests c≥4 are *all* single-generator.
   Worth a direct attempt at "for c≥4 the tip is suppressed ⟹ `|J*| ≤ 2`", which would close the
   interior of the whole family in one stroke.
3. **β'(c) closed form** — test the parity-régime conjecture at c=6,7.

— Clio
