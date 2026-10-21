# Eigenvalue Computation

Eigenvalues describe how a system evolves in time. In a discrete-time linear model,

$$
x_{k+1}=Ax_k,
$$

the state $x_k$ might contain temperatures, displacements, or concentrations. If the initial state is an eigenvector $v$, with $Av=\lambda v$, then

$$
x_k=\lambda^k v.
$$

The eigenvector gives a pattern in the state, called a **mode**. The eigenvalue tells us how that mode changes at each step: it decays if $|\lambda|<1$ and grows if $|\lambda|>1$. A negative or nonreal eigenvalue also introduces changes of sign or phase.

For a continuous-time model $x'(t)=Ax(t)$, the same initial state gives

$$
x(t)=e^{\lambda t}v.
$$

Here $\operatorname{Re}\lambda$ determines exponential growth or decay, while $\operatorname{Im}\lambda$ determines the oscillation frequency. Linear models also approximate nonlinear systems near an equilibrium. Computing eigenvalues therefore helps us identify growing disturbances, slowly decaying modes, and oscillations in scientific and engineering models.

We developed these interpretations in [eigendecomposition](eigendecomposition.md). The task in this chapter is to compute the eigenvalues and the directions associated with them.

## Computational Goals

For $A\in\mathbb C^{n\times n}$, an **eigenpair** is a scalar $\lambda$ and a nonzero vector $v$ satisfying $Av=\lambda v$. Some applications require only a few eigenpairs, such as those governing the fastest growth or slowest decay. Others require the entire **spectrum**, meaning all the eigenvalues of $A$.

We use $^{H}$ for the conjugate transpose. For real matrices, it is the ordinary transpose $^{T}$. Vector symbols are plain lowercase letters, as in the preceding chapters.

For the full spectrum, a useful computational target is the **Schur decomposition**,

$$
A=QTQ^H,
\qquad Q^HQ=I,
$$

where $T$ is upper triangular. Its diagonal entries are the eigenvalues of $A$. The columns of $Q$ are Schur vectors and need not be eigenvectors. Unlike an eigendecomposition $A=X\Lambda X^{-1}$, Schur form exists even when $A$ does not have a basis of eigenvectors. The orthonormal basis also avoids introducing an ill-conditioned change of coordinates.

When working with a real matrix in real arithmetic, the target is real Schur form: an upper block triangular matrix with $1\times1$ and $2\times2$ diagonal blocks. The $2\times2$ blocks represent nonreal conjugate pairs. Eigenvectors, when needed, can be recovered from Schur form. This is the approach used by dense nonsymmetric eigensolvers in [LAPACK](https://www.netlib.org/lapack/lug/node50.html).

## Why Use Iteration?

The definition of an eigenvalue suggests finding the roots of the characteristic polynomial,

$$
p_A(z)=\det(zI-A).
$$

This is useful for theory and small examples. For a general matrix, however, explicitly forming the polynomial and then finding its roots can lose accuracy through rounding its coefficients. Matrix algorithms work directly with $A$ instead.

There is also an algebraic obstacle to a general formula like the quadratic formula. Every polynomial root problem can be written as an eigenvalue problem. To see this, consider a monic polynomial of degree $n\ge2$,

$$
p(z)=z^n+a_{n-1}z^{n-1}+\cdots+a_1z+a_0,
$$

and its **companion matrix**

$$
C=
\begin{pmatrix}
0&1&&&\\
&0&1&&\\
&&\ddots&\ddots&\\
&&&0&1\\
-a_0&-a_1&\cdots&-a_{n-2}&-a_{n-1}
\end{pmatrix}.
$$

````{prf:theorem} Eigenvalues of the Companion Matrix
:label: thm:companion_matrix
The eigenvalues of $C$ are exactly the roots of $p$.
````

````{prf:proof}
Suppose $Cu=zu$, with $u\ne0$. The first $n-1$ rows give

$$
u_{j+1}=zu_j,
\qquad j=1,\ldots,n-1.
$$

Thus $u_j=z^{j-1}u_1$. If $u_1=0$, all entries would vanish, so we can scale $u$ to make $u_1=1$. The last row then gives

$$
-a_0-a_1z-\cdots-a_{n-1}z^{n-1}=z^n,
$$

which is $p(z)=0$.

Conversely, if $p(z)=0$, the nonzero vector $u=(1,z,\ldots,z^{n-1})^T$ satisfies every row of $Cu=zu$. Hence $z$ is an eigenvalue of $C$.
````

The **Abel–Ruffini theorem** states that a general polynomial of degree five or higher has no formula using finitely many arithmetic operations and radicals. The companion matrix shows that the same restriction applies to general eigenvalue formulas. This rules out a general analogue of the quadratic formula; it does not rule out exact solutions for special matrices.

Numerical eigensolvers instead generate successive approximations and stop when an accuracy criterion is met. An eigensolver therefore has two separate questions to answer: does the iteration converge, and how much work is needed to reach a useful approximation?

Accuracy requires the same distinction between algorithm and problem that we used for linear systems. A backward stable eigensolver computes the exact spectrum of a nearby matrix, but sensitive eigenvalues may still have large errors. Iterative convergence and sensitivity both matter; neither fast convergence nor a small backward error alone guarantees accurate eigenvalues.

## From the Power Method to QR Iteration

We begin with the [power method](power_method.md). Repeated multiplication by $A$ can amplify one eigenvector component relative to the others. Under suitable assumptions, normalizing after each multiplication produces an approximation to an eigenvector with a largest-magnitude eigenvalue. This connects the algorithm directly to the time-evolution picture above.

To find more directions, [deflation](deflation.md) removes those already found. [Orthogonal iteration](orthogonal_iteration.md) combines these steps: multiply a set of vectors by $A$, then orthonormalize them using QR. Its convergence analysis explains how ratios of eigenvalue magnitudes control the rate at which the desired subspaces emerge.

Expressing these updates as changes of coordinates leads to [QR iteration](qr_iteration.md), which works toward Schur form. [Hessenberg reduction](upper_hessenberg.md) makes each iteration cheaper, and [shifts](qr_iteration_with_shifts.md) accelerate convergence. We then compare the symmetric and nonsymmetric cases. The following chapter develops methods for computing a few eigenpairs of large, sparse matrices without forming the full Schur decomposition.
