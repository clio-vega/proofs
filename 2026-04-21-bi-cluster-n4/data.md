# Staircase n=4 principal-minor data (compiled 2026-04-21)

Hecke conventions: T_i^2 = (q-1)T_i + q.

Staircase element: Pi_q = (1+T_1)(1+T_2)(1+T_1)(1+T_3)(1+T_2)(1+T_1)
  = B_1 B_2 B_3  where B_j = (1+T_j)(1+T_{j-1})...(1+T_1).

Truncations: Pi_q^(k) = B_1 B_2 ... B_k.  For n=4: k in {1,2,3}, Pi_q^(3) = Pi_q.

Matrix entry: M^(k)[u,v] = [T_u] (T_v * Pi_q^(k))  for u,v in D_n.

## D_4 ordering (lexicographic on one-line notation)

|D_4| = 6.  Entries:

| index | permutation |
|-------|-------------|
| 0 | (2, 1, 4, 3) |
| 1 | (3, 1, 4, 2) |
| 2 | (3, 2, 4, 1) |
| 3 | (4, 1, 3, 2) |
| 4 | (4, 2, 3, 1) |
| 5 | (4, 3, 2, 1) |

## Staircase factor orderings

- `Pi_q^(1)` = (1+T_1)  (word = [1])
- `Pi_q^(2)` = (1+T_1) · (1+T_2) · (1+T_1)  (word = [1, 2, 1])
- `Pi_q^(3)` = (1+T_1) · (1+T_2) · (1+T_1) · (1+T_3) · (1+T_2) · (1+T_1)  (word = [1, 2, 1, 3, 2, 1])

## M^(1)[D_4, D_4]  (6x6)

Rows indexed by u (output), columns by v (input).  M[u,v] = [T_u](T_v Pi_q^(k)).

| u \\ v | (2, 1, 4, 3) | (3, 1, 4, 2) | (3, 2, 4, 1) | (4, 1, 3, 2) | (4, 2, 3, 1) | (4, 3, 2, 1) |
|---|---|---|---|---|---|---|
| (2, 1, 4, 3) | q | 0 | 0 | 0 | 0 | 0 |
| (3, 1, 4, 2) | 0 | q | 0 | 0 | 0 | 0 |
| (3, 2, 4, 1) | 0 | 0 | q | 0 | 0 | 0 |
| (4, 1, 3, 2) | 0 | 0 | 0 | q | 0 | 0 |
| (4, 2, 3, 1) | 0 | 0 | 0 | 0 | q | 0 |
| (4, 3, 2, 1) | 0 | 0 | 0 | 0 | 0 | q |

## M^(2)[D_4, D_4]  (6x6)

Rows indexed by u (output), columns by v (input).  M[u,v] = [T_u](T_v Pi_q^(k)).

| u \\ v | (2, 1, 4, 3) | (3, 1, 4, 2) | (3, 2, 4, 1) | (4, 1, 3, 2) | (4, 2, 3, 1) | (4, 3, 2, 1) |
|---|---|---|---|---|---|---|
| (2, 1, 4, 3) | q*(q + 1) | 0 | 0 | 0 | 0 | 0 |
| (3, 1, 4, 2) | 0 | q*(q + 1) | 0 | q**2 | 0 | 0 |
| (3, 2, 4, 1) | 0 | 0 | q*(q + 1) | 0 | q**2 | q**3 |
| (4, 1, 3, 2) | 0 | q | 0 | 2*q**2 | 0 | 0 |
| (4, 2, 3, 1) | 0 | 0 | q | 0 | 2*q**2 | q**3 |
| (4, 3, 2, 1) | 0 | 0 | q | 0 | q**2 | q**2*(q + 1) |

## M^(3)[D_4, D_4]  (6x6)

Rows indexed by u (output), columns by v (input).  M[u,v] = [T_u](T_v Pi_q^(k)).

| u \\ v | (2, 1, 4, 3) | (3, 1, 4, 2) | (3, 2, 4, 1) | (4, 1, 3, 2) | (4, 2, 3, 1) | (4, 3, 2, 1) |
|---|---|---|---|---|---|---|
| (2, 1, 4, 3) | q**2*(q + 1)*(q + 2) | q**3 | q**4 | 2*q**4 | 2*q**5 | q**6 |
| (3, 1, 4, 2) | q**2 | q**2*(q**2 + 4*q + 1) | q**4 | q**3*(4*q + 1) | q**5 | q**5*(q + 1) |
| (3, 2, 4, 1) | q**2 | q**3 | q**2*(q + 1)*(2*q + 1) | q**4 | q**3*(q**2 + 3*q + 1) | q**4*(q**2 + 3*q + 1) |
| (4, 1, 3, 2) | 2*q**2 | q**2*(4*q + 1) | q**4 | q**3*(5*q + 1) | q**5 | q**5*(q + 1) |
| (4, 2, 3, 1) | 2*q**2 | q**3 | q**2*(q**2 + 3*q + 1) | q**4 | q**3*(q**2 + 4*q + 1) | q**4*(q**2 + 3*q + 1) |
| (4, 3, 2, 1) | q**2 | q**2*(q + 1) | q**2*(q**2 + 3*q + 1) | q**3*(q + 1) | q**3*(q**2 + 3*q + 1) | q**4*(q + 1)*(q + 2) |

## d_k(4) = det M^(k)[D_4]

| k | det (expanded) | det (factored) |
|---|----------------|----------------|
| 1 | `q**6` | `q**6` |
| 2 | `2*q**13 + 7*q**12 + 9*q**11 + 5*q**10 + q**9` | `q**9*(q + 1)**3*(2*q + 1)` |
| 3 | `4*q**26 + 28*q**25 + 85*q**24 + 146*q**23 + 155*q**22 + 104*q**21 + 43*q**20 + 10*q**19 + q**18` | `q**18*(q + 1)**6*(2*q + 1)**2` |

## First-pair identity  d_3 * d_1 ?= q^6 * d_2^2

- LHS (expanded): `q**24*(q + 1)**6*(2*q + 1)**2`
- RHS (expanded): `q**24*(q + 1)**6*(2*q + 1)**2`
- **LHS == RHS?**  `True`

## Eigenvalue decomposition of M^(k)[D_4]

### k = 1

- Characteristic polynomial in t: `(-lambda + q)**6`
- Eigenvalues:
    - q (mult 6)

### k = 2

- Characteristic polynomial in t: `(-lambda + q**2)**3*(-lambda + q**2 + q)*(-lambda + 2*q**2 + q)*(-lambda + q**3 + 2*q**2 + q)`
- Eigenvalues:
    - q*(q + 1) (mult 1)
    - q**2 (mult 3)
    - q*(2*q + 1) (mult 1)
    - q*(q + 1)**2 (mult 1)

### k = 3

- Characteristic polynomial in t: `(-lambda + q**4)*(-lambda**5 + lambda**4*q**6 + 4*lambda**4*q**5 + 14*lambda**4*q**4 + 12*lambda**4*q**3 + 4*lambda**4*q**2 - 7*lambda**3*q**10 - 36*lambda**3*q**9 - 95*lambda**3*q**8 - 122*lambda**3*q**7 - 85*lambda**3*q**6 - 32*lambda**3*q**5 - 5*lambda**3*q**4 + 15*lambda**2*q**14 + 91*lambda**2*q**13 + 256*lambda**2*q**12 + 402*lambda**2*q**11 + 385*lambda**2*q**10 + 235*lambda**2*q**9 + 90*lambda**2*q**8 + 20*lambda**2*q**7 + 2*lambda**2*q**6 - 13*lambda*q**18 - 87*lambda*q**17 - 259*lambda*q**16 - 437*lambda*q**15 - 454*lambda*q**14 - 299*lambda*q**13 - 123*lambda*q**12 - 29*lambda*q**11 - 3*lambda*q**10 + 4*q**22 + 28*q**21 + 85*q**20 + 146*q**19 + 155*q**18 + 104*q**17 + 43*q**16 + 10*q**15 + q**14)`
- Eigenvalue computation failed: It is not always possible to express the eigenvalues of a matrix of size 5x5 or higher in radicals. We have CRootOf, but domains other than the rationals are not currently supported. If there are no symbols in the matrix, it should still be possible to compute numeric approximations of the eigenvalues using M.evalf().eigenvals() or M.charpoly().nroots().

## Known closed forms (from proof 2026-04-18-first-q-shifted-pair.tex)

Expected:
- d_1(4) = q^6
- d_2(4) = q^9 (q+1)^3 (2q+1)
- d_3(4) = q^18 (q+1)^6 (2q+1)^2  =  (d_2(4))^2

| k | matches closed form? |
|---|----------------------|
| 1 | True |
| 2 | True |
| 3 | True |

- n=4 squaring identity  d_3(4) == d_2(4)^2 ?  **True**
  (The corollary of locality: det M^(3) = (det M^(2))^2 at n=4.)

