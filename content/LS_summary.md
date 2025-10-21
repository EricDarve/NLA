# Summary of LS Solution Methods

Choosing the right algorithm for a linear least-squares problem depends on the properties of the matrix $A$, specifically its size, condition number, and rank.

```{list-table} Comparison of Methods for Solving $\text{argmin}_x \|Ax - b\|_2$
:header-rows: 1
:name: lls-methods-summary

* - Method
    - Core Equation(s) to Solve
    - Stability & Conditioning
    - Handles Rank Deficiency?
    - Computational Cost
* - **Normal Equations**
    - 1. Form $C = A^T A$, $d = A^T b$
        2. Solve $Cx = d$
    - **Poor.** $\kappa(A^T A) = \kappa(A)^2$. Numerically unstable for ill-conditioned problems.
    - **No.** Requires $A$ to be full column rank so that $A^T A$ is invertible.
    - $O(mn^2)$. Dominated by forming $A^TA$. Typically the fastest method.
* - **QR Factorization**
    - 1. Factor $A = Q_1R_1$ (Thin QR)
        2. Solve $R_1x = Q_1^T b$
    - **Good.** $\kappa(R_1) = \kappa(A)$. The numerically stable, general-purpose method.
    - **No.** Requires $A$ to be full column rank so that $R_1$ is invertible.
    - $O(mn^2)$. Dominated by the factorization. Slower than Normal Eq.
* - **SVD**
    - 1. Factor $A = U \Sigma V^T$
        2. Solve $x = V \Sigma^\dagger U^T b$
    - **Excellent.** Most numerically robust method, not affected by $\kappa(A)$.
    - **Yes.** This is its key advantage. Finds the unique min-norm solution.
    - $O(mn^2)$. Dominated by the SVD. Generally the slowest of the three.
```