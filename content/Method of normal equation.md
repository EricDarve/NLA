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

[[Least-squares problems]], [[Symmetric Positive Definite Matrices]], [[Cholesky factorization]], [[Conditioning of a linear system]], [[Stability of the Cholesky factorization]]