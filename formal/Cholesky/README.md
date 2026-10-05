# Lean checks for the optional Cholesky material

This project formalizes parts of the argument in
[`content/cholesky.md`](../../content/cholesky.md). It includes a complete proof
of the real inequality behind the **1.1 completion-bound corollary**, together
with matrix and arithmetic lemmas for the counterexample family.

**This is a partial formalization of numerical breakdown. It does not yet prove
that a floating-point Cholesky implementation fails on the proposed family.**
The missing connections are listed below. Compiling this project verifies the
stated Lean theorems, not every claim in the chapter.

## Build and check

The project pins Lean and mathlib to **v4.24.0**. With Lean's `elan` toolchain
manager installed, run from this directory:

```sh
lake exe cache get
lake build
lake env lean Audit.lean
```

See the [mathlib installation instructions](https://leanprover-community.github.io/install/project.html)
for installing the tools. `lean-toolchain` selects Lean, and
`lake-manifest.json` pins the dependency commits. Downloaded dependencies and
compiled files are kept in the ignored `.lake` directory.

`lake build` checks the proofs with warnings treated as errors. `Audit.lean`
prints the axioms used by every theorem in this project. These are limited to
Lean's standard `propext`, `Classical.choice`, and `Quot.sound`. There are no
`sorry` placeholders, added axioms, or proofs by `native_decide` in this project.

## What is checked

| File | Result proved in Lean |
| --- | --- |
| [Completion.lean](Cholesky/Completion.lean) | The exact finite-product inequality behind the constant 1.1, including the square-root sum estimate and logarithmic inequalities. |
| [Schur.lean](Cholesky/Schur.lean) | Positive definiteness of the exact Schur complement; the block-matrix implication used to establish SPD; the all-ones eigenvector of the proposed rounded Schur complement; failure of positive semidefiniteness when $k>9$. |
| [Perturbation.lean](Cholesky/Perturbation.lean) | A negative direction cannot be removed by an error satisfying the stated norm bound. In particular, no Gram factor can satisfy that bound. |
| [Bounds.lean](Cholesky/Bounds.lean) | The scalar estimates giving 57, 8239, the binary64 bound below 532, the lower bound $1156/3$, and the rounding-error margin. |
| [Parameters.lean](Cholesky/Parameters.lean) | At each integer precision $p\ge20$, the largest admissible $q$ exists and satisfies $p-3q\in\{8,10,12\}$. |

### The 1.1 corollary

`Cholesky.meinguet_condition` proves that, for

$$
 n\ge576,\qquad 0<u\le2^{-17},\qquad \kappa\ge1,
 \qquad \frac{11}{10}n\sqrt n\,u\kappa\le1,
$$

we have

$$
 \kappa\left[
 \prod_{j=1}^{n-1}
 \left(1+\frac{u}{1-5u}(\sqrt j+2)\right)-1
 \right]<1.
$$

The Lean expression `n * Real.sqrt n` denotes $n^{3/2}$. The product uses
indices `j = 0, ..., n - 2` and the factor `sqrt (j + 1)`.

This theorem is about real numbers. Meinguet's theorem connecting this product
condition to successful floating-point Cholesky completion remains an external
mathematical result; it is not introduced as an axiom in the Lean project.
Likewise, the variable `κ` is a real number, not a formally defined matrix
condition number.

### The matrix statements

`exactSchur_posDef` proves positive definiteness of

$$
 S=k\eta I_{2k}-\frac{\eta}{8}J_{2k}
$$

for $k>0$ and $\eta>0$, using the quadratic form and Cauchy–Schwarz.
`liftedMatrix_posDef` proves that, whenever $S$ is positive definite,

$$
 \begin{bmatrix}4I&2T\\2T^T&T^TT+S\end{bmatrix}
$$

is positive definite. Its proof expands the quadratic form as
$\|2x+Ty\|_2^2+y^TSy$.

`roundedSchur_ones` proves that the explicit matrix named `roundedSchur`
maps the all-ones vector to $-(k-9)\eta/8$ times itself.
`roundedSchur_not_posSemidef` uses this direction to prove that the matrix is
not positive semidefinite when $k>9$. The name `roundedSchur` refers to the
formula in the chapter. The definition does not execute rounding operations.

`no_gram_factor` proves the perturbation contradiction in a real Hilbert space.
Its hypotheses include a unit vector with quadratic form at most
$-\|A\|/18$ and the backward-error estimate

$$
 \|E\|\le9\frac{s}{1-s}(\|A\|+\|E\|),
 \qquad 0\le s\le\frac{129}{131072}.
$$

The conclusion is that $A+E$ cannot be a Gram operator $R^*R$.
The backward-error estimate is an explicit hypothesis here; it has not been
proved for the floating-point algorithm in this project.

## What remains to prove for the full counterexample

1. Define binary rounding and the precise Cholesky update order. Prove exact
   representability of the input and show that the leading updates produce
   `roundedSchur`, including the binade and exact-subtraction arguments.
2. Formalize the Walsh–Hadamard construction and its Gram identity. This is
   needed to identify the chapter's stored matrix with `liftedMatrix T S`.
3. Derive the sign pattern and componentwise backward-error bound from
   hypothetical successful completion. Prove the operator-norm estimates that
   supply the hypotheses of `no_gram_factor`.
4. Prove the spectral and inverse bounds for the actual matrix, connecting its
   condition number to `conditionUpper x`. The verified inequalities for that
   scalar expression alone are not condition-number theorems.
5. Connect the precision parameters to $x=k^{3/2}u$, and prove that diagonal
   padding preserves both the condition number and numerical failure.

None of these missing claims is silently assumed as a project axiom. They are
separate proof obligations before the full numerical counterexample can be
called formally verified.
