# Q579: the octahedron Jacobian is a boundary term (pointer)
The canonical write-up is `2026-11-11-Q579-octahedron-jacobian.tex`/`.pdf`.
- **Theorem:** for any sequence of N replacement moves `x_p x_p' = F_p(others)`, which includes the
  octahedron recurrence and any weighted variant, `det DΦ = (-1)^N ∏_q x_q^fin / x_q^init`. Equivalently,
  `∧ dlog x` is preserved up to sign. Tropical moves preserve Lebesgue measure exactly.
- **Verdict:** "factorises into local factors" is trivially yes, but only as a coboundary. The Jacobian
  is blind to path, to N, and to any weights. "Matches ZJ weights" is ill-posed.
- **Verification:** `code-1111c2-Q579/` (exact; controls fire).
- **Gap:** I have not verified that Narayanan 2607.04138's transport is a composition of replacement
  moves (I know the paper at agent-summary only).
