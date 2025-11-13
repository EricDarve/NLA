# MINRES

We have seen that GMRES is a general-purpose method for any square matrix, but its cost grows with each iteration. By contrast, the Conjugate Gradient (CG) method is incredibly efficient, but it applies *only* to symmetric positive definite (SPD) matrices.

This raises a natural question: what if we have a matrix $A$ that is **symmetric, but not positive definite (indefinite)?**

* CG is not applicable, as the $A$-norm is not a valid norm and the algorithm is not guaranteed to converge.
* GMRES *is* applicable, but it fails to exploit the symmetry of $A$. We would be paying the full $O(kn)$ cost per iteration, which is unnecessary.

This is the problem solved by the **Minimal Residual (MINRES)** algorithm, developed by Paige and Saunders. MINRES is, in essence, a specialized and highly efficient version of GMRES that leverages symmetry.

## From Arnoldi to Lanczos

The core distinction between GMRES and MINRES is analogous to the difference between the Arnoldi and Lanczos processes:

* **GMRES** is built on the **Arnoldi process** ($A Q_k = Q_{k+1} \underline{H}_k$) for general matrices, producing a Hessenberg matrix $\underline{H}_k$.
* **MINRES** is built on the **Lanczos process** for symmetric matrices. When $A$ is symmetric, the Arnoldi process simplifies *exactly* to the Lanczos process. The resulting Hessenberg matrix $\underline{H}_k$ becomes a tridiagonal matrix, which we will call $\underline{T}_k$.

The Lanczos relation gives us:

$$
A Q_k = Q_{k+1} \underline{T}_k
$$

where $\underline{T}_k$ is a $(k+1) \times k$ *tridiagonal* matrix.

## The MINRES Least-Squares Problem

Like GMRES, MINRES seeks to minimize the 2-norm of the residual at each step:

$$
\min_{y \in \mathbb{R}^k} \|r^{(k)}\|_2 = \min_y \| b - A Q_k y \|_2
$$

We follow the exact same derivation as in GMRES, substituting the Lanczos relation instead of the Arnoldi relation:

$$
\begin{align}
\| b - A Q_k y \|_2
& = \| (\beta e_1) Q_{k+1} - (Q_{k+1} \underline{T}_k) y \|_2 \\[.5em]
& = \left\| Q_{k+1} \left( \beta e_1 - \underline{T}_k y \right) \right\|_2 \\[.5em]
& = \| \beta e_1 - \underline{T}_k y \|_2
\end{align}
$$

Thus, the GMRES problem is reduced from minimizing over a Hessenberg matrix $\underline{H}_k$ to minimizing over a tridiagonal matrix $\underline{T}_k$.

## Efficiency Gains

This reduction in structure is the source of all the efficiency gains. We solve this small least-squares problem just as we do in GMRES: by maintaining a QR factorization using Givens rotations.

$$
G_k^T \cdots G_1^T \underline{T}_k =
\begin{pmatrix}
R_k \\ 0
\end{pmatrix}
$$

However, because $\underline{T}_k$ is sparse (tridiagonal), the Givens rotations required to update this factorization are also very simple. The resulting upper-triangular matrix $R_k$ is also sparse (it has only two upper diagonals).

More importantly, the underlying Lanczos process is a **short-term recurrence**, just like in CG. We no longer need to store the entire $Q_k$ basis or perform the $O(kn)$ orthogonalization.

As a result, the cost per iteration of MINRES is only **$O(\text{nnz} + n)$**.

This is a remarkable result. We have developed an algorithm for general symmetric matrices that has the **exact same asymptotic cost as CG**, but with the broader applicability of GMRES (for symmetric systems).

## Convergence

The convergence of MINRES is governed by the eigenvalue distribution, and for the SPD case, it satisfies a bound similar to that of CG:

$$
\|r^{(k)}\|_2 \le 2 \, \| r^{(0)} \|_2 \left( \frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1} \right)^k
$$

For symmetric indefinite problems, MINRES is the method of choice.

## CG vs. MINRES

* **Conjugate Gradients (CG)** finds the approximation $x^{(k)}_{CG} \in \mathcal{K}_k$ that minimizes the **$A$-norm of the error**:

    $$
    \min_{x \in \mathcal{K}_k} \| x - x^{(k)} \|_A
    $$

* **Minimal Residual (MINRES)** finds the approximation $x^{(k)}_{MINRES} \in \mathcal{K}_k$ that minimizes the **2-norm of the residual**:

    $$
    \min_{x \in \mathcal{K}_k} \| b - A x^{(k)} \|_2
    $$

Even though both are for symmetric systems, both are built on the Lanczos process, and both have the same efficient $O(\text{nnz} + n)$ cost, their iterates ($x^{(k)}_{CG}$ vs. $x^{(k)}_{MINRES}$) are different at each step.

CG guarantees a monotonically decreasing $A$-norm of the error, while MINRES guarantees a monotonically decreasing 2-norm of the residual.

This clarifies the "job" of each algorithm:

* **CG:** The method of choice for **SPD** systems.
* **MINRES:** The method of choice for **symmetric indefinite** systems (where the $A$-norm isn't a norm and CG fails).