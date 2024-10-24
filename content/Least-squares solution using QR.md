Recall our equation from [[Least-squares problems|least-squares]]:
$$
Ax - b \; \perp \; \text{span}(A)
$$
Let's use the [[QR factorization]]: $A = QR$.

Because $R$ is upper triangular we have
$$
\text{span}(A) = \text{span}(Q)
$$
Thus, we get the equivalent linear system
$$
Q^T(Ax - b) = 0
$$
Let's solve for $x$:
$$
Q^T A x = Q^T QR x = Rx = Q^T b
$$
using $A=QR$.

The final solution is
$$
x = R^{-1} Q^T b
$$
Recall the [[Method of normal equation|normal equation]]: 
$$
x = (A^TA)^{-1} A^T b
$$
The condition number is reduced from 
$$
\kappa(A^T A) = \kappa(A)^2
$$
to
$$
\kappa(Q^T A) = \kappa(A)
$$
**This is a huge improvement!**

Note that this method requires $R$ to be non-singular. This is equivalent to saying that $A$ should be full column rank.

The computational cost is $O(mn^2)$.

**Summary**

1. **Orthogonality simplifies the problem**: The orthonormality of $Q$ means that projecting onto $\text{span}(Q)$ (or $\text{span}(A)$) is straightforward. It turns complex geometric relationships into simple algebraic ones.
2. **Upper triangular structure**: Solving systems with upper triangular matrices is simpler because you can solve for one variable at a time, starting from the last equation and moving upwards (back-substitution).
3. **Avoiding squaring the condition number**: By not forming $A^T A$, we prevent the worsening of the condition number, leading to more accurate and stable solutions.
4. **Geometric projection without distortion**: QR factorization allows us to project $b$ onto $\text{span}(A)$ without distorting the space (since $Q$ preserves lengths and angles due to its orthogonality).

[[Least-squares problems]], [[QR factorization]], [[QR using Householder transformations]], [[QR using Givens transformations]], [[Method of normal equation]]