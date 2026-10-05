# Q345 — Lenart–Sottile vs Samuel, 2026-10-05

Two deliberately disjoint code paths, so that agreement between them is evidence.

## Instruments

- `schub.py` — **instrument A**. Schubert polynomials by divided differences from
  `S_{w0}`; structure constants by applying `d_j` along a descent path. Contains
  **no chains, no Monk rule, no Chevalley formula**. Copied (not symlinked —
  a heredoc through a symlink clobbered a script on 10-03) from `code-1004c4-q345/`.
- `chains.py` — **side A**, Samuel's `T_{w/u}` expanded in the root basis, giving
  `N_{w/u}(p)`. From `code-1004c4-q345/`.
- `ls_side.py` — **side B**, NEW. `I_alpha(u,w)` built from the increasing-chain
  definition of `math/0202090` alone. **Never calls `monk_matrix` or `T_tensor`.**
  The only code it shares with instrument A is the definition of `S_n`.

## Scripts, in the order they should be read

| script | what it establishes | log |
|---|---|---|
| `verify_ls.py` | validates the LS enumerator BEFORE any comparison: Prop 3 vs instrument A, Cor 4 vs the structure constants | `verify_ls.out` |
| `compare.py` | the two-sided comparison; symmetry of `T`; `I <= N(p_alpha)`; refutes the content-sum guess | `compare.out` |
| `chi.py` | **the main theorem**: `I_alpha = sum_beta chi_alpha(beta) N(p_beta)`, chi universal | `chi.out` |
| `transition.py` | the converse `eta = N . I^{-1}`, computed from `u=e` only and applied to all `(u,w)` | `transition.out` |
| `pieri_check.py` | the load-bearing identities `S_u h_alpha = sum I_alpha S_w` and `S_u Y^beta = sum N S_w` against instrument A | `n5-pieri.out` |
| `control.py`, `control2.py` | planted negative controls on the LS enumerator | `control2.out` |

## Read the controls

`control.py` planted three and only one fired. `control2.py` is the round-2
diagnosis: three of six fire, and the three silences are each explained — one of
them **confirms** LS's own remark that `b = u(j)` works as well as `b = u(i)`, and
two are vacuous because two *consecutive* labels can never be equal. My first
explanation of that was wrong (I claimed no label can repeat; 598 of 4436 chains
at n=4 do repeat one, non-consecutively). See §7.1 of the paper.

Pure Python, exact integer and `Fraction` arithmetic. No Sage in this container.
