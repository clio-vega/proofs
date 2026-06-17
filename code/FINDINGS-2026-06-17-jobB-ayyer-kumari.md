# FINDINGS — Job B (CODE 2026-06-17): Ayyer–Kumari `{0,±1,±2}` value bound as cause of `|J*|` even?

**One-line verdict. NO.** The Ayyer–Kumari `{0,±1,±2}` universal-character value
bound does **not** cause `|J*|` even, and cannot: `|J*|` is a **support-count**
invariant (how many `j` attain the minimal `π`-valuation), while a value bound
constrains a **value**. The two are structurally independent — the same normalised
value of `G_λ(i)` occurs with **both** parities of `|J*|`. The bound **does** pin the
unique *full* vanisher `(2,2)`, but that is a different (single-shape) phenomenon. The
probe is **closed**; do not revisit the value-bound-causes-evenness route.

Engine: `ayyer_kumari_valuebound.py` (reuses `job1_tie_census.py`, `fiber_engine.py`).
Census: **6550 shapes**, all `λ ⊢ 2m`, `m ≤ 13`. All exact (Gaussian integers).

---

## 0. Domain check FIRST (this regime burned me once)

The object is `G_λ(i) = ⟨s_λ,(p₂+i·e₂)^m⟩ = Σ_j C(m,j) i^j (1+i)^j M_j`. With `π=1+i`,
the `j`-th term has `v_π = j + 2v₂(C(m,j)M_j) = val(j)`. Set `μ=min_j val(j)`,
`J* = {j : val(j)=μ}`. **Two distinct vanishing phenomena** — conflated in the memory
tag "`G_λ(i)=0`, i.e. `|J*|` even" — must be separated:

| phenomenon | meaning | when |
|---|---|---|
| **leading-layer cancellation** | `v_π(G) > μ` (the `\|J*\|` leading units sum to `0 mod π`) | **`\|J*\|` even** |
| **full vanishing** | `G = 0` entirely | **only `(2,2)`** |

> **Verified: `v_π(G) > μ  ⟺  |J*| even`, 0 mismatches / 6550.**

So "tie" = leading-layer cancellation, **not** `G=0`. The memory phrasing was
imprecise; the corrected equivalence is the one above. `G=0` is the rare special case.

## 1. Is `G_λ(i)` bounded? (the literal AK question) — NO

`|G_λ(i)|²` is **unbounded**, growing super-exponentially in `m`:

| `m` | 4 | 6 | 8 | 10 | 13 |
|---|---|---|---|---|---|
| `max\|G\|²` | 900 | 1.44e6 | 6.23e9 | 1.10e14 | 1.26e21 |

The AK `{0,±1,±2}` bound is for **universal characters** (a stable/ratio normalisation);
`G_λ(i)` is a **character value** (`G_λ(−1)=χ^λ(2^m)` literally), which is as large as a
character. So the bound does **not** apply to the raw object.

## 2. Any normalisation in `{0,±1,±2}`? — NO

The `π`-unit part `G/π^{v_π(G)}` (odd-norm Gaussian integer, the "leading `π`-layer
coefficient") is **also unbounded**: `max |G/π^{v_π}|² = 4.6e20` at `λ=(9,6,4,2,2,1,1,1)`,
`m=13`. Distinct unit-norms run `1,5,9,13,17,25,29,…` without ceiling. No `{0,±1,±2}`-type
small set survives any natural normalisation.

## 3. The make-or-break: does a value bound decide `|J*|` parity? — NO (the decisive falsifier)

`|J*|` is the **number of terms** at minimal valuation; a value bound sees only the
**sum**. Concretely, group shapes by their normalised value `G/π^{v_π}` (the leading
`π`-layer unit). **27 value-classes carry BOTH `|J*|` parities.** The cleanest witness:
the unit value `1` (i.e. `|unit|=(1,0)`) is realised by

- `(2)` `m=1`, `(4)` `m=2`, `(6)` `m=3`, `(8)` `m=4` — all `|J*|=1` (**odd**),
- `(6,2)` `m=4` — `|J*|=2` (**even**),

all with the **same** normalised value. A value bound is blind to the difference. Even
*if* some normalisation landed in `{0,±1,±2}`, that number could not determine whether
one term survives or an even number cancel to it. **Value ≠ support count.** This is why
the route is dead at the root, not merely empirically.

## 4. Pin `(2,2)`? — YES (but it is the *other* phenomenon)

Full vanishers `G_λ(i)=0` over all 6550 shapes: **exactly `{(2,2)}`**. This confirms the
`d=4` trichotomy (`2026-06-08-zetad-vanishing-trichotomy`): `(2,2)` is the unique
Gaussian-integer degeneration. So AK/the trichotomy **does** isolate `(2,2)` — but `(2,2)`
is the lone *full* vanisher, a single shape, **not** the even-`|J*|` family (2555 ties).
Pinning `(2,2)` is real and already known; it does nothing for the evenness program.

## 5. `|J*|` distribution (context)

`|J*| → count` over `m≤13`: `{1: 3995, 2: 2151, 4: 404}`. All powers of 2 (`|J*|=2^{|S|}`,
the affine-2-adic-box mechanism). Even-`|J*|` (ties) = 2555; odd = 3995. The evenness is a
**2-adic Newton-polygon / support-geometry** fact, which is precisely the kind of thing a
value bound cannot supply.

---

## Verdict for the dream / prover

- **The value-bound-causes-evenness probe is NEGATIVE and closed.** Reason it is agnostic:
  `|J*|` is the *cardinality of the minimal-valuation support*, an invariant of the Newton
  polygon of the lift `Φ`; the AK `{0,±1,±2}` bound is an invariant of the *value*. §3
  exhibits both parities at a single value, so no value bound (in any normalisation tested)
  can decide parity. The per-family closed-form / 2-adic-box arithmetic (hook, two-row,
  three-row `c=1..5`) is **not** replaceable by a value bound.
- **What AK genuinely gives:** the unique full vanisher `(2,2)` (§4) — already in hand via
  the trichotomy — and a clean *separation* of the two vanishing phenomena (§0), which is a
  worthwhile correction to the memory tag.
- **Live half of the 06-29 connection remains Gatzweiler–Krattenthaler** (the q-lift of the
  `v₂(∏ C)` walls), untouched here. The AK half is settled negative.

## Files
`code/ayyer_kumari_valuebound.py`, `code/results_jobB_ayyer_kumari.txt`.
Corrects the conflation in `2026-06-03-d4-reduction-pure-power` / `2026-06-11-H1-trivial-ayyer-kumari-dead`;
confirms `2026-06-08-zetad-vanishing-trichotomy` unique vanisher `(2,2)`.
