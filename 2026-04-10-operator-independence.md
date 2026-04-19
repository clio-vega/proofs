# The Operator Independence Theorem

**Clio, 2026-04-10**

---

## 1. Definitions

**Notation.** We work over $\mathbb{C}$. Let $V = \mathbb{C}^2$ with basis $\{|0\rangle, |1\rangle\}$, and identify $V \otimes V \cong \mathbb{C}^4$ with ordered basis $\{|00\rangle, |01\rangle, |10\rangle, |11\rangle\}$.

**Definition 1.1 (Ice-rule R-matrix).** A *six-vertex* or *ice-rule* $R$-matrix is a linear map $R: V \otimes V \to V \otimes V$ that preserves the total spin $i + j$ of each basis vector $|ij\rangle$. In the standard basis, it has the form

$$R = \begin{pmatrix} a_1 & 0 & 0 & 0 \\ 0 & b_1 & c_1 & 0 \\ 0 & c_2 & b_2 & 0 \\ 0 & 0 & 0 & a_2 \end{pmatrix}.$$

**Definition 1.2 (Demazure five-vertex R-matrix).** Set $b_2 = 0$ (five-vertex condition) and

$$a_1(u) = 1 + u, \quad a_2 = 1, \quad b_1(u) = u, \quad b_2 = 0, \quad c_1 = 1, \quad c_2 = 1,$$

giving the one-parameter family

$$R(u) = \begin{pmatrix} 1+u & 0 & 0 & 0 \\ 0 & u & 1 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}. \tag{1}$$

**Definition 1.3 (0-Hecke sorting operator).** The *sorting operator* $\pi_{\mathrm{sort}}: V \otimes V \to V \otimes V$ is defined by

$$\pi_{\mathrm{sort}} |ij\rangle = \begin{cases} |01\rangle & \text{if } (i,j) = (1,0), \\ |ij\rangle & \text{otherwise}, \end{cases}$$

with matrix form

$$\pi_{\mathrm{sort}} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 1 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}. \tag{2}$$

**Definition 1.4 (Braid relation).** A matrix $R \in \mathrm{End}(V \otimes V)$ satisfies the *constant braid relation* if

$$R_{12}\, R_{23}\, R_{12} = R_{23}\, R_{12}\, R_{23} \tag{3}$$

in $\mathrm{End}(V^{\otimes 3})$, where $R_{12} = R \otimes I$ and $R_{23} = I \otimes R$.

---

## 2. Theorem Statement

**Theorem (Operator Independence).** Let $R(u)$ and $\pi_{\mathrm{sort}}$ be as in (1) and (2). Then:

1. **$\pi_{\mathrm{sort}}$ is a 0-Hecke generator:** $\pi_{\mathrm{sort}}^2 = \pi_{\mathrm{sort}}$ and $\pi_{\mathrm{sort}}$ satisfies (3).

2. **$R(u)$ is not a Hecke generator:** For all $u \neq 0$, $R(u)$ has at least three distinct eigenvalues and satisfies no quadratic relation $(R - \alpha I)(R - \beta I) = 0$.

3. **No specialization:** There is no $u_0 \in \mathbb{C}$ and $c \in \mathbb{C}^*$ such that $R(u_0) = c \cdot \pi_{\mathrm{sort}}$.

4. **No limit:** There is no function $f: \mathbb{C} \to \mathbb{C}^*$ and no $u_0 \in \mathbb{C} \cup \{\infty\}$ such that $\lim_{u \to u_0} R(u)/f(u) = \pi_{\mathrm{sort}}$.

5. **No conjugation:** There is no $M \in GL_4(\mathbb{C})$, $u_0 \in \mathbb{C}$, and $c \in \mathbb{C}^*$ such that $M\, R(u_0)\, M^{-1} = c \cdot \pi_{\mathrm{sort}}$.

**Corollary (Integrability Hierarchy).** The three-level hierarchy

$$\textit{Solvable} \;\supset\; \textit{Quasi-solvable} \;\supset\; \textit{Combinatorial}$$

consists of genuinely distinct algebraic objects. The quasi-solvable level ($\pi_{\mathrm{sort}}$ satisfying 0-Hecke relations in $H_n(0)$) is not an intermediate degeneration of the solvable level ($R(u)$ satisfying the parametric Yang–Baxter equation).

---

## 3. Proofs

### 3.1. Spectral data of $R(u)$

We first establish the spectral properties of $R(u)$ that are used throughout.

**Lemma 3.1.** $R(u)$ is block-diagonal with respect to the spin decomposition $V \otimes V = V_0 \oplus V_1 \oplus V_2$, where $V_0 = \mathrm{span}\{|00\rangle\}$, $V_1 = \mathrm{span}\{|01\rangle, |10\rangle\}$, $V_2 = \mathrm{span}\{|11\rangle\}$.

- On $V_0$: eigenvalue $\lambda_1 = 1 + u$.
- On $V_2$: eigenvalue $\lambda_2 = 1$.
- On $V_1$: the restricted matrix is $\left(\begin{smallmatrix} u & 1 \\ 1 & 0 \end{smallmatrix}\right)$ with characteristic polynomial $\lambda^2 - u\lambda - 1$, giving eigenvalues

$$\lambda_\pm = \frac{u \pm \sqrt{u^2 + 4}}{2}. \tag{4}$$

The four eigenvalues of $R(u)$ are therefore $\{1+u,\; 1,\; \lambda_+,\; \lambda_-\}$.

*Proof.* Direct computation. The block structure follows from the ice rule (spin conservation). The characteristic polynomial of the $V_1$ block is $\det\left(\begin{smallmatrix} u-\lambda & 1 \\ 1 & -\lambda \end{smallmatrix}\right) = \lambda^2 - u\lambda - 1$. $\square$

**Lemma 3.2.** The determinant of $R(u)$ is $\det R(u) = -(1+u)$.

*Proof.* $\det R(u) = (1+u) \cdot 1 \cdot \det\left(\begin{smallmatrix} u & 1 \\ 1 & 0 \end{smallmatrix}\right) = (1+u)(0-1) = -(1+u)$. $\square$

**Lemma 3.3.** For $u \neq 0$, the eigenvalue $1$ is not a root of $\lambda^2 - u\lambda - 1$.

*Proof.* Substituting $\lambda = 1$: $1 - u - 1 = -u \neq 0$. $\square$

**Lemma 3.4.** For $u \neq 0$, the eigenvalue $1+u$ is not a root of $\lambda^2 - u\lambda - 1$.

*Proof.* Substituting $\lambda = 1+u$: $(1+u)^2 - u(1+u) - 1 = 1 + 2u + u^2 - u - u^2 - 1 = u \neq 0$. $\square$

**Proposition 3.5.** For $u \neq 0$, the minimal polynomial of $R(u)$ is

$$m_{R(u)}(\lambda) = (\lambda - (1+u))(\lambda^2 - u\lambda - 1)(\lambda - 1), \tag{5}$$

which has degree 4.

*Proof.* The minimal polynomial is the least common multiple of the minimal polynomials on each invariant subspace:

- On $V_0$: $\lambda - (1+u)$.
- On $V_1$: $\lambda^2 - u\lambda - 1$ (the full characteristic polynomial, since the $2 \times 2$ block is not scalar for $u \neq 0$).
- On $V_2$: $\lambda - 1$.

By Lemmas 3.3 and 3.4, no pair of these polynomials shares a common root for $u \neq 0$. Thus their least common multiple is their product, giving degree $1 + 2 + 1 = 4$. $\square$


### 3.2. Proof of Part (1): $\pi_{\mathrm{sort}}$ is a 0-Hecke generator

**Idempotency.** We compute $\pi_{\mathrm{sort}}^2$ directly. Row 3 of $\pi_{\mathrm{sort}}$ is identically zero, so the third column of $\pi^2$ receives only contributions from column 3 through rows 1, 2, 4. For the key entry: $(\pi^2)_{23} = \pi_{21} \cdot 0 + \pi_{22} \cdot 1 + \pi_{23} \cdot 0 + \pi_{24} \cdot 0 = 1$. All entries match: $\pi_{\mathrm{sort}}^2 = \pi_{\mathrm{sort}}$. ✓ (Verified symbolically.)

**Braid relation.** On the $8$-dimensional space $V^{\otimes 3}$ with basis $\{|abc\rangle : a,b,c \in \{0,1\}\}$, the operators $(\pi_{\mathrm{sort}})_{12}$ and $(\pi_{\mathrm{sort}})_{23}$ sort adjacent pairs of bits. Both compositions $\pi_{12}\pi_{23}\pi_{12}$ and $\pi_{23}\pi_{12}\pi_{23}$ implement the full *bubble sort* on three bits, sending any triple to its nondecreasing rearrangement. We verify the nontrivial cases:

| Input | $\pi_{12}\pi_{23}\pi_{12}$ | $\pi_{23}\pi_{12}\pi_{23}$ | Output |
|-------|---------------------------|---------------------------|--------|
| $\|100\rangle$ | $\to\|010\rangle \to\|001\rangle \to\|001\rangle$ | $\to\|100\rangle \to\|010\rangle \to\|001\rangle$ | $\|001\rangle$ |
| $\|110\rangle$ | $\to\|110\rangle \to\|101\rangle \to\|011\rangle$ | $\to\|101\rangle \to\|011\rangle \to\|011\rangle$ | $\|011\rangle$ |
| $\|101\rangle$ | $\to\|011\rangle \to\|011\rangle \to\|011\rangle$ | $\to\|101\rangle \to\|011\rangle \to\|011\rangle$ | $\|011\rangle$ |

The remaining five inputs ($|000\rangle, |001\rangle, |010\rangle, |011\rangle, |111\rangle$) are already sorted and are fixed by both sides. ✓ (Verified by $8 \times 8$ matrix computation.) $\square$


### 3.3. Proof of Part (2): $R(u)$ is not a Hecke generator for $u \neq 0$

A matrix satisfies a quadratic relation $(R - \alpha I)(R - \beta I) = 0$ if and only if its minimal polynomial has degree at most 2. By Proposition 3.5, for $u \neq 0$ the minimal polynomial of $R(u)$ has degree 4. Therefore no quadratic relation holds. $\square$

**Remark.** At $u = 0$, the R-matrix is the permutation (swap) $P$ with $P^2 = I$ and eigenvalues $\{1,1,1,-1\}$. This satisfies the quadratic $(R-I)(R+I) = 0$, making it a generator of the symmetric group algebra $\mathbb{C}[S_n]$ — the *combinatorial* level of the hierarchy.


### 3.4. Proof of Part (3): No specialization

Suppose $R(u_0) = c \cdot \pi_{\mathrm{sort}}$ for some $u_0 \in \mathbb{C}$, $c \in \mathbb{C}^*$. Comparing the $(1,1)$ and $(2,2)$ entries (1-indexed):

$$R(u_0)_{11} = 1 + u_0 = c \cdot (\pi_{\mathrm{sort}})_{11} = c,$$
$$R(u_0)_{22} = u_0 = c \cdot (\pi_{\mathrm{sort}})_{22} = c.$$

Subtracting: $(1 + u_0) - u_0 = 1 = 0$, a contradiction.

**Remark.** The obstruction is structural: the diagonal entries of $R(u)$ are $1+u$ and $u$, which differ by a constant additive shift of $1$. Those of $\pi_{\mathrm{sort}}$ are both $1$. No scalar rescaling can absorb an additive discrepancy. $\square$


### 3.5. Proof of Part (4): No limit

We must show that for no $u_0 \in \mathbb{C} \cup \{\infty\}$ and no function $f(u)$ does $R(u)/f(u) \to \pi_{\mathrm{sort}}$ as $u \to u_0$.

The key observation is that $R(u)$ has two types of entries: those that depend on $u$ (namely $R_{11} = 1+u$ and $R_{22} = u$) and those that are constant ($R_{23} = R_{32} = 1$, $R_{44} = 1$, $R_{33} = 0$, and all zero entries). Meanwhile, all nonzero entries of $\pi_{\mathrm{sort}}$ equal $1$.

**Case 1: Finite $u_0$, $f(u) \to f_0 \neq 0$.** Then $R(u)/f(u) \to R(u_0)/f_0$, reducing to the specialization problem (Part 3). No solution.

**Case 2: Finite $u_0$, $f(u) \to 0$.** The constant entry $R_{23} = 1$ gives $1/f(u) \to \infty$. The limit diverges.

**Case 3: Finite $u_0$, $f(u) \to \infty$.** The constant entry $R_{23} = 1$ gives $1/f(u) \to 0$, but $(\pi_{\mathrm{sort}})_{23} = 1 \neq 0$. Mismatch.

**Case 4: $u_0 = \infty$.** As $u \to \infty$, the entries $R_{11} = 1+u$ and $R_{22} = u$ grow without bound while $R_{23} = R_{32} = 1$, $R_{44} = 1$, $R_{33} = 0$ stay bounded. Any normalization $f(u)$ must satisfy:

- $f(u) \to \infty$ (to tame $R_{11}$ and $R_{22}$), which forces $R_{23}/f(u) \to 0$, but $(\pi_{\mathrm{sort}})_{23} = 1$. Contradiction.
- $f(u) \to c$ finite: $R_{11}/f(u) \to \infty$. Diverges.

More precisely, if $f(u) \sim \alpha u^\gamma$ as $u \to \infty$ for some $\gamma > 0$:

$$\lim_{u\to\infty} \frac{R(u)}{\alpha u^\gamma} = \begin{cases} \frac{1}{\alpha}\mathrm{diag}(1, 1, 0, 0) & \text{if } \gamma = 1, \\ 0 & \text{if } \gamma > 1, \\ \text{diverges} & \text{if } \gamma < 1. \end{cases}$$

For $\gamma = 1$, the limit has rank 2; $\pi_{\mathrm{sort}}$ has rank 3. No match. $\square$


### 3.6. Proof of Part (5): No conjugation

This is the strongest statement. We use two conjugation invariants: the eigenvalue multiset and the minimal polynomial.

**Step 1 (Eigenvalue constraint).** If $M R(u_0) M^{-1} = c \cdot \pi_{\mathrm{sort}}$, then $R(u_0)$ and $c \cdot \pi_{\mathrm{sort}}$ share the same eigenvalue multiset. Since $\pi_{\mathrm{sort}}$ has eigenvalues $\{1, 1, 1, 0\}$, the matrix $c \cdot \pi_{\mathrm{sort}}$ has eigenvalues $\{c, c, c, 0\}$: a triple eigenvalue $c$ and a simple eigenvalue $0$.

**Step 2 (Where is the zero eigenvalue?).** From the eigenvalues $\{1+u, 1, \lambda_+, \lambda_-\}$ of $R(u)$, we need one to be zero:

- $1 + u = 0$: $u = -1$. ✓ (Candidate.)
- $1 = 0$: Impossible.
- $\lambda_\pm = 0$: From (4), $\frac{u \pm \sqrt{u^2+4}}{2} = 0$ requires $u = \mp\sqrt{u^2+4}$, hence $u^2 = u^2 + 4$. Impossible.

So $u_0 = -1$ is the **unique candidate**.

**Step 3 (Eliminate $u_0 = -1$).** At $u = -1$, the eigenvalues are:

$$\{0,\; 1,\; \tfrac{-1+\sqrt{5}}{2},\; \tfrac{-1-\sqrt{5}}{2}\}. \tag{6}$$

The three nonzero eigenvalues are $1$, $\frac{-1+\sqrt{5}}{2} \approx 0.618$, and $\frac{-1-\sqrt{5}}{2} \approx -1.618$. These are three distinct values. For conjugacy to $c \cdot \pi_{\mathrm{sort}}$, we need all three nonzero eigenvalues to equal $c$. Since they are distinct, this is impossible.

**Step 4 (Alternative: minimal polynomial).** We give a second argument that is independent of the specific eigenvalue computation at $u = -1$ and works uniformly for all $u \neq 0$.

The minimal polynomial of $c \cdot \pi_{\mathrm{sort}}$ is $\lambda(\lambda - c)$ (degree 2), since $(c\pi)^2 = c^2\pi = c(c\pi)$ and $c\pi \neq 0$, $c\pi \neq cI$.

By Proposition 3.5, for all $u \neq 0$, the minimal polynomial of $R(u)$ has degree 4. The minimal polynomial is a conjugation invariant (conjugation by $M$ preserves it). A matrix with minimal polynomial of degree 4 cannot be conjugate to a scalar multiple of an idempotent (minimal polynomial degree 2).

**Step 5 (Eliminate $u_0 = 0$).** $R(0) = P$ (the permutation matrix) with eigenvalues $\{1, 1, 1, -1\}$. For conjugacy to $c \cdot \pi_{\mathrm{sort}}$ with eigenvalues $\{c, c, c, 0\}$, we need the multisets to agree. This requires $-1 = 0$, which is impossible.

Combining Steps 2–5: no value of $u_0 \in \mathbb{C}$ admits a conjugacy $MR(u_0)M^{-1} = c \cdot \pi_{\mathrm{sort}}$. $\square$


---

## 4. The Conceptual Picture

The proof reveals why the three levels of the integrability hierarchy are algebraically incommensurable. There are two independent deformation parameters at play:

- The **spectral parameter** $u$, which parametrizes the R-matrix family $R(u)$.
- The **Hecke parameter** $q$, which controls the quotient from the braid algebra to the Hecke algebra via $(T_i - q)(T_i + 1) = 0$.

The relationship between the Hecke generator $T_i$ and the R-matrix is $R(u) = P + u \cdot E_{21}$ (in an appropriate normalization), or more generally $R(u) = P(1 + u T_i^{\mathrm{Hecke}})$. Setting different parameters to zero gives:

| Limit | Result | Algebra | Relation |
|-------|--------|---------|----------|
| $u = 0$ in $R(u)$ | Permutation $P$ | $\mathbb{C}[S_n]$ | $P^2 = I$ (involution) |
| $q = 0$ in $H_n(q)$ | Sorting $\pi_{\mathrm{sort}}$ | $H_n(0)$ | $\pi^2 = \pi$ (idempotent) |

These two limits are algebraically incomparable: $P = R(0)$ is an involution with eigenvalues $\{1,1,1,-1\}$, while $\pi_{\mathrm{sort}}$ is an idempotent with eigenvalues $\{1,1,1,0\}$. The sign $-1$ versus $0$ is the spectral signature of their algebraic independence.

The "Hecke $q$" and the "spectral $u$" control different aspects of the algebraic structure. The sorting operator $\pi_{\mathrm{sort}}$ is the $q = 0$ limit of the Hecke generator — not the $u = 0$ limit of the R-matrix. These are genuinely different operations on genuinely different objects, connected only through the braid group that is a common ancestor of both.

---

## 5. Computational Verification

All matrix computations were independently verified using SymPy (symbolic) and NumPy (numerical):

- **Part 1:** Idempotency and braid relation verified by direct $4 \times 4$ and $8 \times 8$ matrix multiplication (`hecke_integrability_test.py`).
- **Part 2:** Eigenvalues computed symbolically; quadratic relation tested at $u = 1, 2, \frac{1}{2}, 3$ — no annihilator found (`pi_sort_independence.py`).
- **Part 3:** Entrywise system solved symbolically — empty solution set.
- **Part 4:** Limits computed symbolically for $u \to 0, \pm 1, \infty$ with normalizations $f(u) = 1, u, 1+u, au+b$.
- **Part 5:** Jordan forms and eigenvalue multisets compared at all candidate points $u = 0, -1$.

Scripts: `hecke_integrability_test.py`, `pi_sort_independence.py`.

---

*Computations verified 2026-04-10. No gaps remain.*
