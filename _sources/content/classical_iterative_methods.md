# Classical Iterative Methods

Splitting methods form the basis of a fundamental class of **stationary iterative methods** used to solve linear systems $A\mathbf{x} = \mathbf{b}$. A stationary iterative method generates a sequence of iterates $\mathbf{x}^{(k)}$ converging to the solution $\mathbf{x} = A^{-1}\mathbf{b}$, where the calculation of $\mathbf{x}^{(k+1)}$ from $\mathbf{x}^{(k)}$ is computationally inexpensive.

A splitting method begins by decomposing the coefficient matrix $A$ into the difference of two matrices, $M$ and $N$ (or sometimes $C$ and $R$, or $M$ and $K$):

$$
A = M - N
$$

Here, $M$ must be chosen such that it is nonsingular and a linear system involving $M$ is easily solvable.

The original equation $A\mathbf{x} = \mathbf{b}$ can then be rewritten as $M\mathbf{x} = N\mathbf{x} + \mathbf{b}$. This leads directly to the core iterative step:

$$
M \mathbf{x}^{(k+1)} = N \mathbf{x}^{(k)} + \mathbf{b}
$$

This can alternatively be expressed as:

$$
\mathbf{x}^{(k+1)} = M^{-1} N \mathbf{x}^{(k)} + M^{-1} \mathbf{b} = R \mathbf{x}^{(k)} + \mathbf{c}.
$$

The convergence of these iterative schemes is determined solely by the **spectral radius** of the iteration matrix $R = M^{-1}N$. These iterative schemes converge for any initial guess $\mathbf{x}^{(0)}$ and $\mathbf{b}$ if and only if $\rho(M^{-1} N) < 1$.

Classical splitting methods are defined by specific choices for $M$:

* **Jacobi Method:** $M = D$, where $D$ is the diagonal of $A$.
* **Gauss-Seidel Method:** $M = D - L$, where $L$ is the strictly lower triangular part of $A$.
* **Successive Over-Relaxation (SOR):** $M = \frac{1}{\omega}D - L$, where $\omega$ is the relaxation parameter. SOR includes Gauss-Seidel as the special case when $\omega=1$.

More advanced techniques, like **Chebyshev iteration**, utilize properties of the spectrum of the iteration matrix to accelerate convergence, positioning them as an effective refinement of the basic stationary approach.

## Motivation: The Problem with Direct Methods for Large Sparse Systems

Splitting methods and other iterative techniques are essential for solving **large sparse linear systems** ($A\mathbf{x} = \mathbf{b}$) because classical **direct solution methods** suffer from catastrophic computational and memory disadvantages in this context.

## The Problem with LU and Cholesky Factorization

Direct methods, such as Gaussian elimination leading to LU factorization (or Cholesky factorization for symmetric positive definite matrices), are highly efficient for dense matrices, typically requiring $O(n^3)$ floating point operations (flops). LU factorization requires $O(n^3)$ flops, and Cholesky factorization requires $n^3/3$ flops.

However, when $A$ is a **large sparse matrix** (often arising from the discretization of Partial Differential Equations, where the number of non-zero entries may be linear, $O(n)$), these direct methods face two critical issues:

1. **High Computational Cost:** While sparsity helps reduce the $O(n^3)$ cost, for truly large systems (especially arising from 3D models), the computational expense remains too high. For a band matrix with semi-bandwidth $w$, the factorization cost is approximately $N w^2/2$ flops, where $N$ is the order of the matrix.
2. **Fill-in:** During the factorization process (e.g., $A=LU$ or $A=LL^T$), new non-zero entries (**fill-in**) are often introduced in positions that were originally zero in $A$. This fill-in drastically increases the memory required to store the factors $L$ and $U$ and complicates the data structures needed to manage the sparse matrix efficiently. In the context of large-scale problems, the demand for computer storage can become unmanageable.

Because of these limitations, iterative methods are practically mandatory for solving systems arising from modern, complex models, such as three-dimensional PDEs.

## Advantages and Limitations of Splitting Methods

The splitting approach provides a viable alternative by reformulating the solution process into an iterative sequence:

### Advantages

*   **Exploitation of Sparsity:** Splitting methods naturally exploit the sparsity of $A$. Since the primary operation required per iteration is a matrix-vector product ($N\mathbf{x}^{(k)}$) and solving a system with $M$ ($M\mathbf{x}^{(k+1)} = \dots$), and $M$ is chosen to be simple (like diagonal or triangular) and $N$ may be sparser than $A$, the computational work per step is minimal.
*   **Computational Simplicity and Parallelism:** These methods are relatively easy to program and require minimal storage (low storage requirement is a main appeal of iterative methods). The updates in methods like Jacobi iteration can often be performed simultaneously or in parallel, making them attractive for parallel computing environments.
*   **Flexibility (Preconditioning):** The matrix $M$ is commonly referred to as the **preconditioning matrix**. The splitting methods themselves, especially in their symmetric variants (SSOR), serve as essential building blocks for more robust and powerful modern techniques, notably the **preconditioned Conjugate Gradient (CG)** method or other Krylov subspace methods.

### Limitations

*   **Convergence Rate and Guarantee:** The method may converge very slowly or, potentially, not converge at all ($\rho(M^{-1}N) \ge 1$). For example, the Jacobi method is guaranteed to converge only if the matrix is strictly row or column diagonally dominant.
*   **Suboptimal Performance:** For sufficiently large problems ($n$ large), basic splitting methods like Jacobi and SOR are generally inferior in performance to highly advanced methods like Fast Fourier Transform (FFT), block cyclic reduction, and Multigrid.
*   **Polynomial Acceleration:** Certain splittings, such as SOR, can lead to iteration matrices where polynomial acceleration (like Chebyshev acceleration) is difficult or yields only minor improvements, especially if the spectrum of the iteration matrix is located on a circle.

## Use Cases, Applications, and Motivation

Splitting methods are highly motivated by applications that lead to massive, sparse linear systems:

*   **Discretization of PDEs:** A primary application is solving systems arising from finite difference or finite element discretizations of Partial Differential Equations, particularly elliptic types. In these cases, the coefficient matrix $A$ is typically large and highly sparse, with a well-defined structure (e.g., matrices consisting of a few nonzero diagonals).
*   **Preconditioning:** The most significant modern use case for classical splitting methods is as preconditioners. The matrix $M$ is sought to be a simple, easily invertible approximation of $A$ such that the preconditioned system $M^{-1} A \mathbf{x} = M^{-1} \mathbf{b}$ has a much smaller condition number $\kappa(M^{-1}A)$ than $\kappa(A)$. A small condition number leads to fast convergence when combined with Krylov subspace acceleration methods (like Conjugate Gradient or GMRES).
*   **Building Blocks for Multilevel Methods:** Basic iterative methods remain relevant as they often serve as "smoothers" within sophisticated multilevel methods like Multigrid, where their simplicity and ability to quickly damp high-frequency error components are leveraged.