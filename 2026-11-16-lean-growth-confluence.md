# LEAN 2026-11-16: growth-diagram confluence (lem:D, half (ii))

**Target.** Half (ii) of Lemma `lem:D` (inward = outward switching) in
`proofs/2026-11-13c2-Q585-G4-evacuation.tex`. It is the abstract confluence argument: any two monotone flip
sequences from the initial labelled path to the final one give the same labels, and both equal the recursion-defined η.

**Project.** `lean/tworow_d4_kernel`, module `TworowD4Kernel/GrowthConfluence.lean`. It is imported by the root, so it is
in `defaultTargets`. Repo: github.com/clio-vega/tworow-d4-kernel.

## Encoding
- The grid is `[0,m]×[0,p]`. A path runs from (0,p) to (m,0) and is stored as a word in `{R,D}` with label function `f : ℕ → V`.
  Position t of the word is the paper's η_t, so the encoding is literally the paper's (word, chain) pair.
- Initial path: `D^p R^m`. Final path: `R^m D^p`.
- `Flip L`: a flip rewrites `a ++ D::R::b` to `a ++ R::D::b` and sets `f(|a|+1) := L f(|a|) f(|a|+1) f(|a|+2)`. The new vertex is the
  far (upper-right) corner of the square.
- `eta L p g`: `η(0,j)=g(p−j)`, `η(i,0)=g(p+i)`, `η(i,j)=L η(i−1,j) η(i−1,j−1) η(i,j−1)`. It is defined by recursion on i+j.
- Invariant `Good`: the word has p letters `D`, every vertex t carries `eta (vtx t)`, and every non-vertex position carries `g`.
  It is about the current path only. The claim "each vertex is created exactly once" is **not** used.

## Sorry-free results (0 sorries in the module)
| declaration | content |
|---|---|
| `good_of_reachable`, `reachable_label_eq_eta` | along any flip sequence, each vertex on the path carries η |
| `confluence` | two sequences ending on the same word give **equal** label functions |
| `confluence_final` | the special case used in Lemma D (end at `R^m D^p`) |
| `good_unflip` | backward moves preserve the invariant if `L a (L a b c) c = b` |
| `flipAt_sound` | the executable flip used for the lesions is a `Flip` |

## Lesions (kernel `decide`)
1. **`lesion_step_indexed`** (load-bearing). The rule `L s a b c = s` reads the step number. On the 2×2 grid, rows-first
   `[1,0,2,1]` gives labels `[0,1,3,2,0]` and columns-first `[1,2,0,1]` gives `[0,2,3,1,0]`. Both runs are legal and both end at `RRDD`.
   Positive control: `control_two_orders_agree` with the step-blind rule `a+2b+3c+1`.
2. **`lesion_unflip`**. A backward move with the non-involutive rule `L = b+1` on the 1×1 grid: "flip" gives `[0,1,0]`, while
   "flip, unflip, flip" gives `[0,3,0]`. The positive half is `good_unflip`: with an involutive rule, backward moves are harmless.
3. **`lesion_initial`** (control). Different initial labels give different final labels.

## Axioms (one run, sorryAx control in the same file)
```
confluence              : [propext, Classical.choice, Quot.sound]
confluence_final        : [propext, Classical.choice, Quot.sound]
reachable_label_eq_eta  : [propext, Classical.choice, Quot.sound]
good_unflip, flipAt_sound : [propext, Quot.sound]
lesion_*, control_*     : no axioms
sorry_control           : [sorryAx]
```

## Build and test
- `lake build`: green, 3493 jobs.
- `lake test`: green (~10 min, almost all of it the pre-existing checks). Three new `#guard`s shadow the control and lesions 1 and 2.
- Planted fault (ordB→ordA in the lesion-1 guard), checked standalone (`lake env lean`): the real guard passes and the planted one fails (`did not evaluate to true`, rc=1). A full planted `lake test` timed out at 10 min and was killed; the test file was restored and diffed against the passing version.
- Lean commit: tworow-d4-kernel `134f9fb` (upstream HEAD verified).

## Not formalised
- **Half (i)**: a forward/reverse jdt slide is a product of elementary switches. This is tableau-specific. The registry parent
  `H1-G4-rotation-map-is-switching` therefore stays at `proved`. The new child `G4-growth-confluence-formal` is `lean-verified`.
- The general-m,p rows-first and columns-first sequences are not constructed as explicit flip lists. `confluence_final` covers
  every pair of flip sequences, these two included, and the 2×2 instance is checked concretely.
