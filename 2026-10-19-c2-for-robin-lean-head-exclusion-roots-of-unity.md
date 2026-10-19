# For Robin — LEAN 2026-10-19 c2: `cor:domain` formalised, and the characteristic hypothesis is spent in only one place

**Lean file:**
<https://github.com/clio-vega/tworow-d4-kernel/blob/main/TworowD4Kernel/HeadExclusionRootsOfUnity.lean>
**Session note (full lesion register, build evidence, declaration inventory):**
<https://github.com/clio-vega/proofs/blob/main/2026-10-19-c2-lean-head-exclusion-roots-of-unity.md>
**This note:**
<https://github.com/clio-vega/proofs/blob/main/2026-10-19-c2-for-robin-lean-head-exclusion-roots-of-unity.md>
**Paper being formalised:**
<https://github.com/clio-vega/proofs/blob/main/2026-10-19-typeB-gauge-theorem.tex> (`cor:domain`, L727–L752)

Sorry-free. 36 named declarations (6 defs, 30 theorems), 0 sorries, axioms exactly
`[propext, Classical.choice, Quot.sound]`.

## What it says

Yesterday's prove session found that the **type-`B` head factor is `(1−q)(v−1)`** where
type `C`'s is `(1−q)(v²−1)`. Since a `τ`-fixed point has `v² = 1`, type `C` fires at
both `±1` and type `B` only at `+1` — which is why type-`B` head exclusion **fails at
`−1`**, and why the theorem's domain is cut out by that type's own root system.

The Lean file defines

```lean
def headFactor (q v : K) (h : ℕ) : K := (1 - q) * (v ^ h - 1)
def tauFiring (q : K) (h : ℕ) : Set K := {v | v ^ 2 = 1 ∧ headFactor q v h = 0}
```

with `h` a **free natural-number parameter**, and proves

```lean
theorem tauFiring_eq_inter (hq : q ≠ 1) :
    tauFiring q h = tauFixedLocus K ∩ boundaryLocus K h
theorem tauFiring_of_odd  (hq : q ≠ 1) (ho : Odd h)  : tauFiring q h = ({1} : Set K)
theorem tauFiring_of_even (hq : q ≠ 1) (he : Even h) : tauFiring q h = tauFixedLocus K
```

Type `B` and type `C` are `h = 1` and `h = 2` of these.

## The thing worth your attention

**The mechanism is the parity of `h` and nothing else.** On the `τ`-fixed locus
`v^h = v^(h % 2)`, so the type-`B`/type-`C` dichotomy is the odd/even dichotomy at its
two smallest representatives.

That matters because of a retraction I made yesterday. My own `2026-10-17-typeB`
`thm:rhoB` glossed `h` as *"the length of the boundary root: `2t_c` in type `C`, `t_c`
in type `B`"*, and that is wrong — `|2e_m|² = 4`, `|e_m|² = 1`, neither is `h`. `h` is
the exponent in the boundary **character** `w_m^h`. I retracted it in prose in
`rem:h`.

**The formalisation is evidence about that claim rather than a restatement of it.** Lean
admits `h` only as a `ℕ` parameter of every declaration, and the proof shows it enters
solely through `h % 2`. A quantity that was secretly a type-dependent *length* could
not have that form. This is the one place in the development where the `h`-uniformity
is real, and the type checker is what certifies it.

## The lesion result, which corrects my own brief

My brief for this session told me: *"work over a field of characteristic ≠ 2, state
that hypothesis, and lesion it — it should be irremovable, and if it is not, I have
formalised something weaker than I think."*

It **is** irremovable — but for exactly **one** of the two statements, and not the one
I was thinking of. The whole characterisation is **characteristic-free**:
`tauFiring_eq_inter`, `fires_at_tauFixed_iff`, `tauFiring_of_odd`, `tauFiring_of_even`,
`tauFiring_typeB`, `tauFiring_typeC`, `tauFixedLocus_eq` all hold over `𝔽₂`, where they
correctly degenerate (`+1 = −1`, the `τ`-fixed locus is one point, both types' loci are
`{1}`).

`(-1) ≠ 1` is spent **only** on the theorem that *separates the two types*:

```lean
theorem tauFiring_typeB_ssubset_typeC (hq : q ≠ 1) (hchar : (-1 : K) ≠ 1) :
    tauFiring q 1 ⊂ tauFiring q (2 : ℕ)
```

and the collapse in characteristic 2 is itself a **theorem**, not a failed
counterexample search:

```lean
theorem tauFiring_typeB_eq_typeC_of_charTwo [CharP K 2] (hq : q ≠ 1) :
    tauFiring q 1 = tauFiring q (2 : ℕ)
```

So: *a genericity hypothesis belongs to the theorem that makes the distinction, not to
the theorem that describes the object.* That is word for word what the 10-17
`CollapseMultiplier` session found about char-2 in the type-C locus lemma — reproduced
on an unrelated statement, one day and one lemma later. I only caught it because
10-18's rule says to lesion a binder in **every** theorem that carries it; from either
declaration alone the verdict is invisible.

Other verdicts, each with a witness in the file rather than a failed rebuild:
`q ≠ 1` irremovable; `v² = 1` irremovable in `fires_at_tauFixed_iff` but **relocated**
out of `fires_iff_pow_eq_one`; `[IsDomain K]` irremovable (over `ZMod 6` with `q = 3`
the head factor vanishes at `v = 2`, where `2² = 4 ≠ 1`, both factors nonzero, product
`12 = 0`); `h ≥ 1` removable; `v ≠ 0` absent and vacuously satisfied. And `Field` is
more than the theorem needs — `CommRing + IsDomain` suffices, nothing is inverted.

## What I did NOT formalise, so the registry parent stays `proved`

1. **`lem:head` itself** — the vanishing of `T^B_σ` — which needs the six-line
   `T`-operator, head selections and `cor:heads`. No declaration mentions `T`, `σ` or
   `r`. What is here is the head **factor**'s zero locus, which is `lem:head`(b)'s
   discriminator.
2. **`def:D`'s root-system input**: that the boundary root of `B_N` is `e_m` and of
   `C_N` is `2e_m`, hence `h_B = 1` and `h_C = 2`. Those enter as the numerals `1` and
   `2`.
3. **`cor:domain`(1)'s root-system half** — that `{w_m = −1}` is not a root hyperplane
   of `B_N`. Only the arithmetic half is formalised (`-1 ∉ {w : w = 1}`). So the "G2
   surviving is not luck" reading still rests on an unformalised root-system fact, and
   I would not want it quoted as Lean-verified.

The `lean-verified` grade went on a new registry **child**,
`lean-head-exclusion-roots-of-unity`, under `adjunction-domain-is-the-root-system`,
which stays `proved`. Registry: `proofs/registry/motegi-typeBC-multiplication.json`,
validated clean (`OK: ... is valid`).

## One instrument fault worth passing on

A background build reported **exit code 0 on a build that failed.** I had wrapped it as
`lake build > log 2>&1; echo "LAKE_EXIT=$?" >> log`, so the shell's status was the
`echo`'s, not Lake's. The log said `LAKE_EXIT=1` and `error: build failed`. I only
caught it because I grep the log rather than trusting the status — but the healthy and
broken cases produce the *same* status, which is the same shape as the stalled-vs-broken
`lake build` problem from 10-18. If you ever see me report a green build from a wrapper's
exit code, distrust it.

Also: `lake build -j1` does not exist in Lake 4.30 (`unknown short option '-j'`, and
there is no `--jobs` in `lake --help`), so "serialise the build" has no flag behind it —
it just means don't run two builds at once.
