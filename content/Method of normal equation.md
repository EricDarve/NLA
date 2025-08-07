The method of normal equation consists in [[Least-squares problems|solving]]
$$
(A^TA) \, x = A^T b
$$
The solution is:
$$
x = (A^TA)^{-1} A^T b
$$
The matrix $A^TA$ is [[Symmetric Positive Definite Matrices|SPD]]. So the system can be solved using [[Cholesky factorization|Cholesky]].

This method is best for very tall skinny $A$.

One of the main drawbacks is that the [[Conditioning of a linear system|condition number]] grows very quickly! Indeed we can prove that
$$
\kappa(A^T A) = \| A^T A\|_2 \: \|(A^TA)^{-1}\|_2 = \kappa(A)^2
$$
So the condition number grows much faster than $\kappa(A)$.

This method requires $A^TA$ to be non-singular. This is equivalent to saying that $A$ should be full column rank.

The computational cost is $O(mn^2)$.

**Intuitive explanation**
- $A^T A$: This product represents the "correlation" of A's columns with each other. It captures how the columns of A interact and overlap.
- $A^T b$: This term represents the "correlation" of A's columns with the target vector b. It tells us how much each column of A contributes to explaining b.
- $(A^T A)^{-1}$: Inverting $A^T A$ is like "decorrelating" the columns of A. It accounts for any redundancy or overlap in A's columns.
- Final multiplication: $(A^T A)^{-1} A^T b$ combines the decorrelated version of A with its correlation to b, giving us the optimal coefficients x.

[[Least-squares problems]], [[Symmetric Positive Definite Matrices]], [[Cholesky factorization]], [[Conditioning of a linear system]], [[Stability of the Cholesky factorization]]