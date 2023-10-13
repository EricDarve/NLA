The [[Cholesky factorization|Cholesky]] factorization algorithm is very stable. The entries of $L$ cannot grow. **This is true even if we do not [[Row pivoting|pivot]].**

This can be seen from the equations:
$$
A = L L^T, \qquad \sum_j l_{ij}^2 = a_{ii}
$$
Therefore the entries of $L$ are [[Backward error analysis for LU|bounded]] by the square root of the diagonal entries of $A.$

[[Stability of the LU factorization]], [[Backward error analysis for LU]], [[Row pivoting]], [[Cholesky factorization]], [[Existence of the Cholesky factorization]]

