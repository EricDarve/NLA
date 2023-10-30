We compare the computational cost of performing the QR iteration on unsymmetric vs [[Hermitian and symmetric matrices|symmetric matrices]]. As we expect, the symmetric case is significantly faster.

### Unsymmetric

- Cost of [[Upper Hessenberg form for the QR iteration|upper Hessenberg form:]] $O(n^3)$
- [[QR iteration for upper Hessenberg matrices|QR iteration step:]] cost $O(n^2)$ using Givens transformations.

### Symmetric

- Cost of [[Upper Hessenberg form for the QR iteration|upper Hessenberg form:]] $O(n^3)$

From $Q^T A Q = H$, we get that $H$ is [[Hermitian and symmetric matrices|symmetric.]] So $H$ is Hessenberg and symmetric. So it is tri-diagonal symmetric.

Consider now the steps in the [[QR iteration]]. Recall that 
$$
T_{k+1} = U_{k+1}^H T_k U_{k+1}.
$$
So $T_k$ remains [[Hermitian and symmetric matrices|symmetric.]] It is also upper Hessenberg, so it remains tri-diagonal symmetric.

- The cost of each [[QR iteration for upper Hessenberg matrices|QR iteration step]] is just $O(n)$. 
- The first step is expensive, but after that, each iteration is extremely cheap.