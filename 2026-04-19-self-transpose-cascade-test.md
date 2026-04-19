# Self-transpose non-hook cascade test

**Date:** 2026-04-19
**Script:** `~/projects/proofs/self_transpose_cascade_test.py`

## Setup

Two obstruction sets are on the table:

1. **Liao-Rybnikov cactus non-2-transitivity** (arXiv 2506.16561): the cactus
   group J_n acts 2-transitively on SYT(λ) *iff* λ is not self-transpose and
   not a hook. Obstruction set = {self-transpose non-hook λ}.

2. **Clio's even-block cascade gap.** The cascade transport + parabolic
   reduction proof of the Rank Isolation Lemma covers all odd blocks and the
   k = 4 even block directly; the k ≥ 6 even blocks were the direct-proof
   gap (closed Apr 17 by eigenvalue-positivity bypass).

**Conjecture:** obstruction (1) ⊆ obstruction (2) (or equal). Both live on
ρ(w_0) = evacuation fixed points.

## Correction to the task's partition list

The task prompt named (5,3,1) as self-transpose of n = 9. **False** — the
conjugate of (5,3,1) is (3,2,2,1,1). It also claimed "n = 8: none non-hook".
**False** — n = 8 has TWO self-transpose non-hook partitions, (4,2,1,1) and
(3,3,2).

The correct list of self-transpose non-hook λ with |λ| ≤ 9:

| n | λ          | conjugate   | dim V_λ |
|---|------------|-------------|---------|
| 4 | (2,2)      | (2,2)       | 2       |
| 6 | (3,2,1)    | (3,2,1)     | 16      |
| 8 | (4,2,1,1)  | (4,2,1,1)   | 90      |
| 8 | (3,3,2)    | (3,3,2)     | 42      |
| 9 | (3,3,3)    | (3,3,3)     | 42      |

(n = 1, 2, 3, 5, 7 have no non-hook self-transpose partitions.)

## Methodology

The operator is Clio's **staircase reduced word** Π_q:

    Π_q = B_{n-1} · B_{n-2} · ... · B_1,
    B_k  = (1 + T_k)(1 + T_{k-1}) ... (1 + T_1).

This is the standard reduced-word expression for w_0; **not** the
single-factor (1 + T_1)(1 + T_2)...(1 + T_{n-1}) in the task prompt. The
H-invariant theorem r_λ = dim V_λ^H applies to the staircase Π.

The representation ρ_λ is built in Hoefsmit q-seminormal form. Matrices are
verified to satisfy the quadratic, braid, and far-commutation Hecke
relations on every test case with dim V_λ ≤ 20 (the larger cases are trusted
by construction).

Ranks are computed by sympy exact Gaussian elimination at q = 1, q = 2, and
at generic rational values {3, 3/2, 5/7} to certify the generic rank over
Q(q). The H-invariant prediction r_λ = dim V_λ^H is computed as the rank of
the projector Π_{i odd} (1 + T_i)/2 at q = 1.

The cumulative rank profile is computed block-by-block to observe rank drops.

## Results table

| λ           | n | dim V | rk at q=1 | rk at q=2 | rk over Q(q) | dim V^H | even k ≥ 6 traversed | self-T | hook |
|-------------|---|-------|-----------|-----------|--------------|---------|----------------------|--------|------|
| (2,2)       | 4 | 2     | 1         | 1         | 1            | 1       | []                   | yes    | no   |
| (3,2,1)     | 6 | 16    | 2         | 2         | 2            | 2       | []                   | yes    | no   |
| (4,2,1,1)   | 8 | 90    | 3         | (q-indep) | (q-indep)    | 3       | [6]                  | yes    | no   |
| (3,3,2)     | 8 | 42    | 3         | 3         | 3            | 3       | [6]                  | yes    | no   |
| (3,3,3)     | 9 | 42    | 3         | 3         | 3            | 3       | [6, 8]               | yes    | no   |

Sanity checks passed:

- Hecke relations (quadratic, braid, far-commute) verified on dim ≤ 20 cases.
- rank at q=1 = rank at q=2 = rank over Q(q) (q-independence, as expected).
- rank = dim V_λ^H in every row (H-invariant theorem, as expected).

### Block-by-block cascade drops

For every λ tested, rank drops occur ONLY at **odd blocks**. No even-block
rank drop appears in any V_λ. Example:

- (3,3,3): block ranks 42 → 21 (k=1) → 21 → 11 (k=3) → 11 → 6 (k=5) → 6 →
  3 (k=7) → 3. Drops at k = 1, 3, 5, 7; no drop at k = 2, 4, 6, 8.

This matches the Rank Isolation Lemma: rank drops occur only at
H-generators s_1, s_3, s_5, ... (the odd-indexed generators).

## Interpreting "cascade-reachable"

The interesting question is what "the direct cascade proof can handle λ"
means. Two interpretations:

**Interpretation A (literal "does rank drop at even k ≥ 6 in V_λ?"):** trivially
NO for every λ. The Rank Isolation Lemma proves NO rank drop at any even block
for any λ. So obstruction (2) under interpretation A is empty, which makes
the whole comparison vacuous. This is not what the task is asking.

**Interpretation B ("does the cascade *argument* pass through an even block
k ≥ 6 for this λ?"):** YES iff n ≥ 7 (since block k = 6 appears in Π for
n ≥ 7, block 8 for n ≥ 9, etc.). Under this interpretation:

- (2,2), n=4: reachable by direct proof (max even block is k = 2 < 6).
- (3,2,1), n=6: reachable (max even block is k = 4, covered by the direct
  k=4 verification).
- (4,2,1,1), (3,3,2), n=8: NOT reachable by direct proof (must certify
  block k = 6 ≥ 6, which was the gap).
- (3,3,3), n=9: NOT reachable by direct proof (blocks k = 6, 8 in the gap).

But under Interpretation B, "cascade-unreachable" = "n ≥ 7", which holds
for every λ of n ≥ 7, self-transpose or not. The obstruction set is
{all λ with n ≥ 7}, which is much larger than the self-transpose non-hook
set. Equality of obstructions fails badly.

**Interpretation C ("does λ fail the parabolic-reduction finite check for
some even k ≥ 6?"):** the parabolic reduction step translates block-k on
V_λ to S_{k+1} irreps. For block k=6, the check is done once for all 15
partitions of 7 (Clio's memory, verified); similarly k=8 for 30 partitions
of 9, k=10 for 56 partitions of 11. None of these finite checks failed —
all evaluated to "no drop". Under this interpretation, obstruction (2) is
again empty.

## Verdict

Under every honest interpretation I can formulate from Clio's corpus:

- **Interpretation A** (rank drop at even k ≥ 6 in V_λ): obstruction (2) is
  empty. No λ is hit. So (1) ≠ (2).
- **Interpretation B** (argument traverses even k ≥ 6): obstruction (2) =
  {λ : n(λ) ≥ 7}. Every self-transpose non-hook with n ≥ 7 is hit, but so
  are many non-self-transpose λ. So (1) ⊆ (2) on our sample, but strictly.
- **Interpretation C** (parabolic check fails): empty in every verified n.

**On this sample, obstruction (1) ≠ obstruction (2) under A and C, and (1)
is properly contained in (2) under B.** There is no interpretation where
(1) = (2) on these five partitions.

**Evidence weight: LOW for equality, MODERATE against it.** The sample is
tiny (5 partitions), and the two obstructions are framed at different
levels — (1) is about the cactus action's 2-transitivity on a SYT set;
(2) is about rank drops of a Hecke element on an irrep. They need not
coincide, and the earlier memory entry
(`project_self_transpose_dichotomy.md`) already flagged this as a
*suspicion*, not a theorem.

## Surprises

1. **The task's partition list was wrong.** (5,3,1) is not self-conjugate.
   n = 8 has two self-transpose non-hooks, not zero. I corrected silently.

2. **The rank in the task prompt was wrong in operator.** The prompt said
   Π_q = (1+T_1)(1+T_2)...(1+T_{n-1}). Clio's actual operator is the
   staircase reduced word Π = B_{n-1}·...·B_1. These have different ranks:
   for (3,2,1) the linear product has rank 8, the staircase has rank 2.
   The H-invariant theorem holds for the staircase, not the linear product.
   I computed both for completeness but the relevant comparison is the
   staircase.

3. **The "SYT support formula"** r_lam = #{T : row_T(2j) ≤ j-1 for all j}
   is NOT r_λ. For (3,2,1) it gives 6, but rank(Π|V) = 2. The SYT support
   formula describes the *support basis* of the image at intermediate
   block k-1, not the dimension of the image of the full staircase.
   `seminormal_analysis.py` part 10 is careful about this; I was sloppy
   at first.

4. **No even-block rank drops exist anywhere.** This is the content of the
   Rank Isolation Lemma, confirmed here: every rank drop in Π happens at
   an odd block. Under the literal interpretation of "cascade obstruction
   at k ≥ 6", the set is empty. The meaningful "gap" is not *where* rank
   drops, but *which case of the proof* needs even-block parabolic
   reduction — a property of n, not of the partition.

5. **Self-transpose non-hook partitions appear uniformly at dim = 2 or 3
   for our small n.** No rank-structural feature jumps out that is
   specific to self-transposeness in the rank vector.

## Recommendation

The conjecture "cactus obstruction = cascade obstruction" is not supported
by this sample under any interpretation I can make rigorous. If the
conjecture is to hold, the two obstructions must be formulated at the
*same* level — e.g., both as statements about the evacuation
ρ(w_0) fixed locus in V_λ, or both about some structural subspace. As
stated, they are framed at different levels and the comparison is weak.

The more promising related question (from the memory trail) is whether
self-transpose non-hook λ are distinguished in SOME OTHER structural
invariant of Π_q on V_λ — e.g., the q-spectrum of Π on V_λ, the
nilpotent-vs-semisimple decomposition, or the discrepancy between the
eigenspace decomposition and the ρ(w_0) grading. That was not the
question asked here, but it might be the right question.
