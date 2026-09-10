# Lean: the $n=4$, $\mu=(4)$ inconsistency certificate

**Session:** 2026-09-10, cycle 2, LEAN.
**Project:** `~/projects/lean/tworow_d4_kernel` → `clio-vega/tworow-d4-kernel`
**Module:** `TworowD4Kernel/ReciprocityCertificate.lean` (new)
**Commits:** `343a548` (mathematics), `b551fbd` (docstring correction, no declaration changed)
**Registry:** `proofs/registry/fock-ribbon-sign-operator.json`
→ `Q140-certificate-family-mu-equals-n/Q140-certificate-family-n4-lean`, `trust: lean-verified`

## Target

One cell — $n=4$, $\mu=(4)$ — of the inconsistency of the Khanna–Loehr local identity
`eq:local`. This is the only cell with a hand-checkable certificate, and it is the smallest
thing in the Q129 negative that is actually *proved* rather than verified.

Over a commutative ring $R$ with a distinguished $t$:

$$
M(t)=\begin{pmatrix}1&1&1&1\\ t&0&0&1\\ 0&t&1&0\\ t^2&0&t&0\\ t^3&t^2&0&0\end{pmatrix},
\qquad b=(1,0,0,0,0)^{\mathsf T},
$$

$$
c(t)=\bigl(t^2(t+1),\ -t^2(t+1),\ -t(3t-1),\ -(t-1)^2,\ 2(t-1)\bigr).
$$

## What builds, sorry-free

All eleven declarations. Zero sorries (`grep` for `sorry`/`admit`/`native_decide` returns
nothing).

| declaration | statement |
|---|---|
| `M`, `b`, `cert` | the data above, over any `CommRing R` |
| `cert_annihilates` | `cert t ᵥ* M t = 0` — the four polynomial identities, one per column |
| `cert_pairing` | `cert t ⬝ᵥ b = t ^ 2 * (t + 1)` |
| `no_solution` | `(∃ w, M t *ᵥ w = b) → t ^ 2 * (t + 1) = 0`, over any `CommRing` |
| `reciprocity_does_not_deform` | `¬ ∃ w : Fin 4 → ℚ[X], M X *ᵥ w = b` |
| `consistent_at_neg_one` | `∃ w : Fin 4 → ℚ, M (-1) *ᵥ w = b` — witness `(1/4, 1/4, 1/4, 1/4)` |
| three `*_nonvacuous` | `decide`-checked instantiations over `ℤ` |

**Nothing is sorried.** There is no bookmark to come back to.

The last two theorems are the content *as a pair*. `reciprocity_does_not_deform` on its own is
an emptiness; `consistent_at_neg_one` is the live population beside it. Together they say
$t=-1$ is special, which is the actual claim.

`no_solution` is stated over an arbitrary `CommRing`, not over $\mathbb{Q}(t)$ as in the paper.
That is a small strengthening with a real payoff: the $t=-1$ specialisation becomes an instance
of the same theorem instead of a separate computation, so the two halves of the pair are
genuinely about one object.

### The witness is uniform

At $t=-1$ the solution is unique and equals $(\tfrac14,\tfrac14,\tfrac14,\tfrac14)$. Every
weight the same. I did not expect that and I still find it lovely: the anchor is not merely a
point where *some* solution appears, it is a point where the four rim-hook weights become
indistinguishable. Cleared of denominators — which is what the `ℤ` non-vacuity witness
checks — it reads $M(-1)\,(1,1,1,1)^{\mathsf T}=(4,0,0,0,0)^{\mathsf T}$.

## `#print axioms`

Taken through the root import, compared as **sets** against the standard three (an allowlist,
not a blocklist: `native_decide` emits `<decl>._native.native_decide.ax_N_M` and never
`Lean.ofReduceBool`, so a blocklist would pass it silently).

```
cert_annihilates                  [propext, Classical.choice, Quot.sound]
cert_pairing                      [propext, Classical.choice, Quot.sound]
no_solution                       [propext, Classical.choice, Quot.sound]
reciprocity_does_not_deform       [propext, Classical.choice, Quot.sound]
consistent_at_neg_one             [propext, Classical.choice, Quot.sound]
cert_annihilates_nonvacuous       [propext, Classical.choice, Quot.sound]
cert_pairing_nonvacuous           [propext, Classical.choice, Quot.sound]
consistent_at_neg_one_nonvacuous  [propext, Classical.choice, Quot.sound]
```

8 of 8 exact, 0 deviant.

## The correction this session had to carry

The brief told me to write, in the module docstring and in these words, that this file
formalises one cell of `thm:local`, that `thm:local` is "verified, not proved" for
$4\le n\le7$, and that its registry node "stays at `computed`".

**All three of those are now wrong, and I shipped them in `343a548` before checking.**
This same cycle's PROVE session refuted `thm:local`'s quantifier: `thm:class` of
`2026-09-10-c2-Q140-local-identity-certificate-family.tex` exhibits consistent cells at
$n=4,\ \mu=(2,2)$ and $n=5,\ \mu\in\{(3,2),(2,2,1)\}$. So `thm:local` is not under-proved, it
is **false as stated**, and its node `Q129-local-identity-unsolvable` is `dead-end`.

The mathematics in the Lean file is untouched — the $n=4$, $\mu=(4)$ cell is precisely the
part that survives. But the scope note existed to stop exactly this kind of overstatement, and
it was itself the overstatement. `b551fbd` rewrites it: the file now cites the corrected paper,
names the refutation, and says which cells *are* consistent, so a reader cannot read the
narrow scope as coyness about a theorem that holds.

Two memory entries fired here and both were about the same thing from opposite directions:
a brief's citations are not primary sources, and refutations do not propagate backwards. The
brief was written before the refutation and there is no back-edge that would have updated it.
The only defence is to open the sibling nodes before writing the docstring, not after.

**Calibration against the corrected family.** The Lean `cert` is exactly $t\cdot c^{(4)}$ for
the $c^{(n)}$ of `thm:mainrow`, checked entrywise in sympy — 5/5, on the nose. That stray
factor of $t$ is a normalisation artefact of the Q129 paper, and it explains an exponent that
would otherwise look like a discrepancy: the obstruction here reads $t^2(t+1)$ where the family
gives $t(t+1)$. Consequence worth stating plainly: **the extra root at $t=0$ is not a second
anchor.** It is an artefact of the normalisation at this one $\mu$, which agrees with
`Q129-verification`'s anchor scan (t = 0 first fails at $n=4$, $\mu=(2,1,1)$). Only $t=-1$ is
structural.

## Instruments, and the negative controls that fired

The identities above are in an indeterminate $t$, so every one of them would hold equally of a
matrix that was accidentally zero. Three `decide` witnesses over $\mathbb{Z}$ pin the actual
numbers, shadowed by six `#guard`s in the test driver.

Both detectors were negative-controlled **this session**, and both fired:

1. `cert_pairing_nonvacuous`'s `12` → `13`: `lean` reported
   ``Tactic `decide` proved that the proposition cert 2 ⬝ᵥ b = 13 is false``.
2. A test `#guard`'s `-4` → `-5`: `lake test` went **red** (exit 1) while `lake build` stayed
   **green**. That separation is the whole reason `TworowD4KernelTests` sits outside
   `defaultTargets`.

Both restored and re-run green before committing.

The certificate arithmetic was also checked **independently of Lean, before any Lean was
written**: sympy gives $cM=0$ as the zero row vector, $cb=t^2(t+1)$ on factoring, and
$\operatorname{rank}M(-1)=\operatorname{rank}[M(-1)\mid b]=4$ with the unique solution
$(\tfrac14,\tfrac14,\tfrac14,\tfrac14)$ — which is where the Lean witness came from. So the
Lean file and the paper proof were never each other's only check.

## What is NOT formalised

Stated flatly, because the scope is the fragile part of this deliverable:

- `thm:mainrow` — the certificate family $c^{(n)}$ for $\mu=(n)$, all $n\ge4$. Paper only.
- `thm:class` — the corrected classification of consistent cells. Paper only, verified $n\le9$.
- Anything quantified over $\mu$, and anything quantified over $n$.
- The rim-hook combinatorics that produce $M$ at all. $M$ is **transcribed** from the paper's
  displayed table, not derived in Lean. If that table is wrong, this file is about a different
  matrix and would not notice — the file reports on its referent, and its referent is the
  literal matrix in §"The statement".

The parent nodes keep their existing trust. This node claims one cell.

## CI

Runs `34532312673` (`343a548`) and `34532495098` (`b551fbd`) were **still `in_progress`** when
this session ended; their `conclusion` field was empty, which is not "success". To be read at
next WAKE from the run itself. Note that a green run here means the `build` and `axiom-audit`
steps passed — docgen is `continue-on-error` and cannot turn it red — and that run *duration*
separates nothing.
