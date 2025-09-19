# The Sherman-Morrison-Woodbury Formula

The **Sherman-Morrison-Woodbury formula** provides an efficient and powerful method for calculating the inverse of a matrix that has been updated by a low-rank matrix. Its primary advantage is that it allows you to find the new inverse without having to perform a full, computationally expensive matrix inversion from scratch, provided you already know the inverse of the original matrix. 💡

## The Sherman-Morrison Formula (Rank-One Update)

The simplest case of this identity is the **Sherman-Morrison formula**, which applies to a **rank-one update**. A rank-one update means we are perturbing an invertible matrix $A$ by adding an outer product of two vectors, $uv^T$.

If $A$ is an invertible matrix and $u, v$ are vectors, then the inverse of the updated matrix $A + uv^T$ is given by:

$$(A + uv^T)^{-1} = A^{-1} - \frac{A^{-1}uv^T A^{-1}}{1 + v^T A^{-1} u}$$

The remarkable efficiency of this formula comes from the denominator, $1 + v^T A^{-1} u$. Since $v^T$, $A^{-1}$, and $u$ are a row vector, a matrix, and a column vector, respectively, their product is a **scalar**. This means we avoid a complex matrix inversion and instead only need to perform a simple scalar division.

## Application: Solving a Linear System

In practice, we often use the formula to efficiently solve a linear system of the form $(A + uv^T)x = b$, especially when $A^{-1}$ is known or systems with $A$ are easy to solve.

To compute $x = (A + uv^T)^{-1} b$, we can apply the formula and group the operations intelligently:

$$x = \left( A^{-1} - \frac{A^{-1}u(v^T A^{-1})}{1 + v^T A^{-1} u} \right) b = A^{-1}b - A^{-1}u \left( \frac{v^T A^{-1} b}{1 + v^T A^{-1} u} \right)$$

This leads to the following efficient computational steps:
1.  Solve for an intermediate vector $y = A^{-1}b$.
2.  Solve for another intermediate vector $z = A^{-1}u$.
3.  Compute the scalar value in the denominator: $\beta = 1 + v^T z$.
4.  Compute the scalar value in the numerator's right part: $\alpha = v^T y$.
5.  Combine these results to find the final solution: $x = y - z (\frac{\alpha}{\beta})$.

This procedure replaces a full $O(n^3)$ matrix inversion with a few matrix-vector multiplications (or system solves, which are typically $O(n^2)$) and vector operations, resulting in significant computational savings.

## The Woodbury Matrix Identity (Generalization)

The Sherman-Morrison formula is a special case of the more general **Woodbury matrix identity**, which handles updates of a higher rank. The formula is stated as:

$$(A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^{-1}U)^{-1}VA^{-1}$$

Here, the update is $UCV$, where $U$ is an $n \times k$ matrix, $C$ is a $k \times k$ matrix, and $V$ is a $k \times n$ matrix. The key advantage is that instead of inverting the large $n \times n$ matrix on the left, we only need to invert the much smaller $k \times k$ matrix $(C^{-1} + VA^{-1}U)$ on the right. This is extremely beneficial when $k$ is much smaller than $n$.