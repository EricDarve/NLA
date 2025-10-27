# QR iteration for upper Hessenberg matrices

The primary advantage of first converting a matrix to upper Hessenberg form is the significant reduction in computational cost for each step of the QR iteration. A single QR iteration for a dense matrix costs $O(n^3)$, but for an upper Hessenberg matrix, this cost is reduced to **$O(n^2)$**.

The process works as follows:

1.  **Step 1: QR Factorization ($H = QR$)**
    * The goal is to transform the upper Hessenberg matrix $H$ into an upper triangular matrix $R$.
    * Because $H$ is already *almost* upper triangular (it only has one subdiagonal), this is done very efficiently using $n-1$ Givens transformations.
    * Each Givens rotation $G_k$ is applied from the left to zero out a single subdiagonal element $H[k+1, k]$.
    * The complete $Q$ matrix (which is orthogonal) is the product of all these individual Givens rotations: $Q^T = G_{n-1} \dots G_2 G_1$. Thus, $R = Q^T H$.
    * Applying each of the $n-1$ rotations to the matrix takes $O(n)$ operations, leading to a total cost of **$O(n^2)$** for this step.

2.  **Step 2: Re-multiplication ($H_{new} = RQ$)**
    * To complete the iteration, the matrices are multiplied in the reverse order: $H_{new} = RQ$.
    * Since $R = Q^T H$, we can write this as $H_{new} = (Q^T H) Q$. However, to maintain efficiency, we don't form $Q$ explicitly.
    * Instead, we apply the saved Givens rotations (specifically, their transposes, $G_k^T$) one by one from the *right* to the $R$ matrix: $H_{new} = R G_1^T G_2^T \dots G_{n-1}^T$.
    * This step also involves $n-1$ rotations, each costing $O(n)$ operations, for a total cost of **$O(n^2)$**.

## Summary

Crucially, this two-step process **preserves the upper Hessenberg form**, so the new matrix $H_{new}$ is also upper Hessenberg.

This allows the next iteration to be performed just as efficiently. While a single iteration is $O(n^2)$, finding all eigenvalues typically requires $O(n)$ iterations, leading to a total time cost of **$O(n^3)$** for the entire algorithm, a significant improvement over the $O(n^4)$ required by the naive QR iteration on a full matrix.