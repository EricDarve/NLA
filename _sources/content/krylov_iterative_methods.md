# Krylov Subspace Methods

## Searching for an Optimal Solution

In the previous chapters, we explored classical iterative methods like Jacobi, Gauss-Seidel, and SOR. These "splitting" methods provide a simple and intuitive way to approach the linear system $Ax=b$, but their convergence can be slow and highly dependent on the properties of the matrix $A$. We now turn our attention to a more powerful and sophisticated family of techniques: **Krylov subspace methods**.

These methods represent a fundamental shift in strategy. Instead of a simple fixed-point iteration, we will build a sequence of approximate solutions $x^{(k)}$ by searching within a carefully constructed, expanding subspace. This subspace, known as the **Krylov subspace**, is the heart of the entire approach.

**Definition: The Krylov Subspace**

Given a matrix $A$ and a starting vector (typically the initial residual $r^{(0)} = b - Ax^{(0)}$), the $k$-th Krylov subspace ${\mathcal K}_k$ is the space spanned by the first $k$ vectors in the sequence $\{r^{(0)}, Ar^{(0)}, A^2r^{(0)}, \ldots\}$:

$$
{\mathcal K}_k(A, r^{(0)}) = \text{span}\{ r^{(0)}, A r^{(0)}, A^2 r^{(0)}, \ldots, A^{k-1} r^{(0)} \}
$$

The core idea of Krylov methods is to find an $x^{(k)}$ in the affine subspace $x^{(0)} + {\mathcal K}_k$ that is "optimal" in some sense. Assuming $x^{(0)} = 0$ for simplicity, our goal is to find $x^{(k)} \in {\mathcal K}_k(A, b)$ that best approximates the true solution $x$.

## The Central Challenge: What is "Optimal"?

This brings us to the central question: **How do we define the "optimal" solution** within the Krylov subspace?

Let's say we have an orthogonal basis $Q_k$ for ${\mathcal K}_k$. Our approximate solution would then be $x^{(k)} = Q_k y_k$ for some vector of coefficients $y_k$.

The most natural and desirable choice would be to find the $x^{(k)}$ that is *closest* to the true, unknown solution $x$. This means we would try to minimize the 2-norm of the **true error** $e^{(k)} = x - x^{(k)}$:

$$
y_k = \underset{y \in \mathbb{R}^k}{\rm argmin} \, \| x - Q_k y \|_2
$$

The solution to this least-squares problem is the orthogonal projection of $x$ onto the subspace, given by $y_k = Q_k^T x$.

Here we hit a wall. To compute $Q_k^T x$, we would need to know $x$. But $x = A^{-1}b$ is the very solution we are trying to find! This "best" possible solution is intractable.


## Two Tractable Paths to Optimality

Instead, Krylov methods define alternative, *computable* optimality criteria that we *can* solve. The choice of criterion depends on the properties of the matrix $A$, and different choices give rise to different methods.

In this chapter, we will focus on the two most important and widely used methods.

1.  **The Conjugate Gradient Method (CG)**
    * **For:** Symmetric Positive Definite (SPD) matrices.
    * **The Idea:** For SPD matrices, we can define a special "energy" norm, known as the **A-norm**: $\|v\|_A = \sqrt{v^T A v}$. It turns out that finding the $x^{(k)} \in {\mathcal K}_k$ that minimizes the **$A$-norm of the true error** is a tractable problem. This insight leads to the celebrated Conjugate Gradient method, a remarkably efficient algorithm that (in exact arithmetic) guarantees convergence in at most $n$ steps.

2.  **The Generalized Minimal Residual Method (GMRES)**
    * **For:** General, nonsingular matrices.
    * **The Idea:** When $A$ is not SPD, the $A$-norm is not available. Instead, we take a different approach. While we cannot measure the error $e^{(k)}$, we *can* measure the **residual** $r^{(k)} = b - Ax^{(k)}$. The GMRES method seeks the $x^{(k)} \in {\mathcal K}_k$ that minimizes the 2-norm of this residual.

        $$
        \text{Find } x^{(k)} \in {\mathcal K}_k(A, b) \quad \text{such that} \quad \| b - A x^{(k)} \|_2 \text{ is minimized.}
        $$
        
    This is a standard linear least-squares problem (in $k$ dimensions) that we can solve at each step.

We will begin by deriving the Conjugate Gradient method, leveraging the special properties of SPD matrices, before moving on to the more general (and more complex) framework of GMRES.