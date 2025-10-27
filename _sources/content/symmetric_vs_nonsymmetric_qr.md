# Symmetric vs Unsymmetric QR Iteration

We discuss the differences between the reduction to Hessenberg form, the shifted QR iteration, and their behavior for symmetric versus unsymmetric real matrices. We also provide additional details on the QR iteration algorithm with shifts for both cases.

The standard, high-performance approach for finding all eigenvalues of a dense matrix is a **two-phase process**:

1.  **Phase 1 (Direct):** Use a *direct* (non-iterative) method to reduce the matrix to a simpler form that is unitarily similar. This step is $O(n^3)$.
2.  **Phase 2 (Iterative):** Use an *iterative* method (the shifted QR iteration) on this simpler form to find the eigenvalues.

The differences between the symmetric and unsymmetric cases affect every part of this process.

## Phase 1: Reduction to Upper Hessenberg Form

The goal of this phase is to introduce as many zeros as possible using a *finite* number of similarity transformations (e.g., Householder reflections) to make the subsequent iterative phase cheaper.

### Unsymmetric Real Matrices

* **Target Form:** **Upper Hessenberg form**.
* **Structure:** A matrix $H$ is upper Hessenberg if all entries below the first sub-diagonal are zero ($H_{ij} = 0$ for $i > j+1$). It's the "closest" you can get to upper triangular in a finite number of steps.

    $$ H = \begin{pmatrix}
    \times & \times & \times & \times \\
    \times & \times & \times & \times \\
    0 & \times & \times & \times \\
    0 & 0 & \times & \times
    \end{pmatrix}
    $$

* **Computational Cost:** This reduction requires $\approx \frac{10}{3}n^3$ floating-point operations (flops).

### Symmetric Real Matrices

* **Target Form:** **Tridiagonal form**.
* **Structure:** If the original matrix $A$ is symmetric ($A = A^T$), any similarity transformation $H = Q^T A Q$ will also be symmetric. A **symmetric upper Hessenberg matrix is, by definition, tridiagonal**.

    $$ T = \begin{pmatrix}
    \times & \times & 0 & 0 \\
    \times & \times & \times & 0 \\
    0 & \times & \times & \times \\
    0 & 0 & \times & \times
    \end{pmatrix}
    $$

* **Computational Cost:** This reduction is significantly cheaper, requiring $\approx \frac{4}{3}n^3$ flops.

This initial reduction is a *huge* win for the symmetric case. The target matrix $T$ has only $O(n)$ non-zero elements, compared to the $O(n^2)$ elements in the unsymmetric Hessenberg matrix $H$.

## Phase 2: The Shifted QR Iteration

This is the iterative phase where we apply QR steps to the reduced matrix ($H$ or $T$) to find the eigenvalues.

### Unsymmetric Real Matrices

* **Eigenvalues:** The eigenvalues can be **real or complex conjugate pairs**.
* **The Problem:** A simple shift like $\mu = H[n,n]$ will fail to converge if the eigenvalue it's "chasing" is part of a complex pair. The shift (which is real) will oscillate and never settle.
* **The Algorithm (Francis QR Step):** The standard solution is the **double-shift QR algorithm**.
    1.  It uses *two* shifts, $\mu_1$ and $\mu_2$, which are chosen as the eigenvalues of the bottom-right $2 \times 2$ submatrix of $H$.
    2.  If the eigenvalues are complex (e.g., $a \pm bi$), this choice perfectly captures them.
    3.  A "bulge chasing" technique performs this double-step implicitly, using only real arithmetic, which is very clever and stable.
* **Cost per Iteration:** This "bulge chase" on a Hessenberg matrix costs **$O(n^2)$** operations.
* **Convergence:** The convergence is generally quadratic with the Francis double-shift: 

    $$
    \epsilon_{k+1} \approx c \, \epsilon_k^2
    $$ 

    for some constant $c < 1$, where $\epsilon_k$ is the error at iteration $k$.
* **Final Form:** The iteration converges to the **Schur form**—an upper triangular (or block-upper-triangular) matrix $S$. The eigenvalues are on the diagonal (or in $2 \times 2$ blocks for complex pairs).

### Symmetric Real Matrices

* **Eigenvalues:** The eigenvalues are **always real**.
* **The Algorithm (Wilkinson Shift):** Because the eigenvalues are real, we do not need a double-shift. We can use a single, real shift.
    1.  While the $\mu = T[n,n]$ shift works, a much more robust and faster-converging choice is the **Wilkinson shift**. It's the eigenvalue of the bottom $2 \times 2$ submatrix that is *closer* to $T[n,n]$.
    2.  This shift provides a safeguard against some pathological cases where the simple $T[n,n]$ shift might fail.
* **Cost per Iteration:** Since the matrix is tridiagonal, applying one step of the QR algorithm (with "bulge chasing") is astonishingly cheap. It only costs **$O(n)$** operations.
* **Convergence:** The convergence for the symmetric case with a Wilkinson shift is **globally and cubically convergent**. This is incredibly fast; it's one of the most powerful results in numerical linear algebra:

    $$
    \epsilon_{k+1} \approx c \, \epsilon_k^3
    $$

    for some constant $c < 1$.
* **Final Form:** The iteration converges to a **diagonal matrix** $D$. (A symmetric upper triangular matrix is, by definition, diagonal). The eigenvalues are the diagonal entries.

## 3. Summary: Cost and Convergence

Here is a side-by-side comparison of finding all eigenvalues.

| Feature | Symmetric Case (e.g., `np.linalg.eigh`) | Unsymmetric Case (e.g., `np.linalg.eig`) |
| :--- | :--- | :--- |
| **Eigenvalues** | Always **real** | Can be **real or complex conjugate pairs** |
| **Phase 1: Form** | **Tridiagonal** ($O(n)$ non-zeros) | **Upper Hessenberg** ($O(n^2)$ non-zeros) |
| **Phase 1: Cost** | $\approx \frac{4}{3}n^3$ flops | $\approx \frac{10}{3}n^3$ flops |
| **Phase 2: Shift** | Single real shift (e.g., Wilkinson) | Double shift (Francis QR step) |
| **Phase 2: Cost / Iteration**| **$O(n)$** | **$O(n^2)$** |
| **Phase 2: Convergence**| **Cubic** (extremely fast) | **Quadratic** |
| **Phase 2: Total Cost** | **$O(n^2)$** | **$O(n^3)$** |
| **Overall Total Cost** | $O(n^3) + O(n^2) = \mathbf{O(n^3)}$ | $O(n^3) + O(n^3) = \mathbf{O(n^3)}$ |
| **Final Form** | **Diagonal** (Schur form is diagonal) | **Upper Triangular** (Schur form) |

**Key Takeaway on Cost:**

While both algorithms have the same asymptotic complexity of $O(n^3)$, the **symmetric case is significantly faster**.
* The $O(n^3)$ reduction step has a smaller constant ($\frac{4}{3}$ vs. $\frac{10}{3}$).
* The $O(n^3)$ iterative phase for the unsymmetric case is completely *avoided* in the symmetric case, which has an almost-negligible $O(n^2)$ iterative cost.
* Therefore, the total cost for the symmetric problem is dominated *only* by the initial reduction, whereas both phases contribute significantly to the cost of the unsymmetric problem.