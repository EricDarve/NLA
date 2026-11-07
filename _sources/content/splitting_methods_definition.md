# Splitting Methods and Convergence Theory

The [chapter introduction](classical_iterative_methods.md) described iteration as a sequence of inexpensive corrections. A matrix splitting makes that idea precise. We first derive the update, then ask when repeated updates converge.

## From a Splitting to an Iteration

Consider $Ax=b$, where $A$ is a real nonsingular $n\times n$ matrix. Write $x$ for the exact solution and $x^{(k)}$ for the approximation after $k$ iterations, starting from a chosen $x^{(0)}$. The superscript counts iterations, not powers.

Choose a nonsingular matrix $M$ whose linear systems are inexpensive to solve, and write

$$
A=M-N.
$$

This is a **splitting** of $A$. Since the exact solution satisfies $Mx=Nx+b$, we define the iteration by

$$
Mx^{(k+1)}=Nx^{(k)}+b.
$$

Each step uses the current approximation on the right and solves a system with $M$ to obtain the next approximation. Equivalently,

$$
\begin{aligned}
x^{(k+1)}&=Rx^{(k)}+c,\\
R&=M^{-1}N=I-M^{-1}A,
\qquad c=M^{-1}b.
\end{aligned}
$$

The matrix $R$ is the **iteration matrix**. The method is **stationary** because $R$ and $c$ are fixed throughout the iteration. In practice, we solve with $M$; these formulas do not require forming $M^{-1}$.

The update also has a useful correction form. Define the residual $r^{(k)}=b-Ax^{(k)}$. Subtracting $Mx^{(k)}$ from the iteration equation gives

$$
M\bigl(x^{(k+1)}-x^{(k)}\bigr)=r^{(k)}.
$$

The exact correction would satisfy $A(x-x^{(k)})=r^{(k)}$. Replacing $A$ by $M$ makes this correction cheaper to compute. Choosing $M=A$ would solve the original problem in one step, but would require the very solve we are trying to avoid. Simpler choices are cheaper, but may need many steps or fail to converge.

## How the Error Evolves

Let $e^{(k)}=x^{(k)}-x$. Subtract the exact equation $Mx=Nx+b$ from the iteration equation:

$$
\begin{aligned}
M\bigl(x^{(k+1)}-x\bigr)&=N\bigl(x^{(k)}-x\bigr),\\
e^{(k+1)}&=Re^{(k)}.
\end{aligned}
$$

Repeated substitution gives

$$
\boxed{e^{(k)}=R^ke^{(0)}.}
$$

Convergence for every initial guess is therefore equivalent to $R^k\to0$. This connects the linear-system iteration to the behavior of matrix powers studied in the eigenvalue chapter.

## The Convergence Criterion

If an error lies along an eigenvector of $R$ with eigenvalue $\lambda$, each step multiplies it by $\lambda$. That error decays exactly when $|\lambda|<1$. This suggests checking the largest eigenvalue magnitude, called the **spectral radius**:

$$
\rho(R)=\max_j|\lambda_j(R)|.
$$

:::{prf:theorem} Convergence of Splitting Methods
:label: thm:splitting_convergence

Let $A=M-N$, with $A$ and $M$ nonsingular. The iteration converges to $A^{-1}b$ for every initial guess and every right-hand side if and only if

$$
\boxed{\rho(M^{-1}N)<1.}
$$

This result does not require the iteration matrix to be diagonalizable.
:::

:::{prf:proof}
We have reduced the question to when $R^k\to0$, where $R=M^{-1}N$.

**Necessity.** For any eigenpair $Rz=\lambda z$ with $z\ne0$,

$$
R^kz=\lambda^kz.
$$

If $R^k\to0$, this vector must tend to zero, which requires $|\lambda|<1$. Hence $\rho(R)<1$.

**Sufficiency when $R$ is diagonalizable.** Write $R=X\Lambda X^{-1}$. If $\rho(R)<1$, every diagonal entry of $\Lambda^k$ tends to zero, so

$$
R^k=X\operatorname{diag}(\lambda_1^k,\ldots,\lambda_n^k)X^{-1}
\longrightarrow0.
$$

This proves the result for diagonalizable matrices. The extension below shows why the same criterion holds without an eigenvector basis.
:::

### Extension to Nondiagonalizable Matrices (Optional)

:::{note}
This argument completes the proof in the general case. It uses a Schur decomposition and a change of scale between coordinates. You can skip it on a first reading and use the convergence criterion as stated.
:::

:::{prf:proof}
Assume $\rho(R)<1$. A diagonalization may not exist, but a complex Schur decomposition always does:

$$
R=QTQ^H.
$$

As in the [Schur decomposition](eigendecomposition.md) studied earlier, $Q$ is unitary, $Q^H$ denotes its conjugate transpose, and $T$ is upper triangular. Its diagonal entries are the eigenvalues of $R$. The entries above the diagonal couple the coordinates; we make those entries small by rescaling.

For $0<\delta<1$, set

$$
D=\operatorname{diag}(1,\delta,\ldots,\delta^{n-1}),
\qquad B=D^{-1}TD.
$$

Then $b_{ij}=t_{ij}\delta^{j-i}$ for $i\le j$. The diagonal is unchanged, while every entry above it tends to zero as $\delta\to0$. Consequently, the maximum absolute row sum satisfies

$$
\|B\|_\infty
=\max_i\left(
|\lambda_i|+\sum_{j>i}|t_{ij}|\delta^{j-i}
\right)
\longrightarrow\rho(R)<1.
$$

Choose and fix a sufficiently small $\delta$ so that $\|B\|_\infty<1$. Repeated use of the [matrix norm inequality](operator_norms.md) $\|BC\|_\infty\le\|B\|_\infty\|C\|_\infty$ gives

$$
\|B^k\|_\infty\le\|B\|_\infty^k\longrightarrow0.
$$

Since $Q$, $D$, and their inverses are fixed,

$$
R^k=QD B^kD^{-1}Q^H\longrightarrow0.
$$

Thus the error tends to zero for every initial guess, including when $R$ is not diagonalizable.
:::

## What the Criterion Tells Us

The eigenvalues determine whether all initial errors eventually decay. They do not guarantee that the error decreases at every step in the usual norm. Coupling between coordinates can cause temporary growth even when $\rho(R)<1$.

A spectral radius close to one can make convergence slow. The choice of splitting therefore matters both for the cost of each step and for the number of steps needed. The next sections examine this tradeoff for [Jacobi](jacobi_method.md), [Gauss-Seidel](gauss_seidel_method.md), and [SOR](sor_method.md).
