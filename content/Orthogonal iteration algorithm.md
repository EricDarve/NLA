Pseudo-algorithm for the [[Orthogonal iteration|orthogonal iteration algorithm:]]

```julia
while not_converged:
    Qk = A * Qk
    Qk, Rk = qr(Qk)
```

Consider:
$$
T_k = Q_k^H A Q_k
$$
$T_k$ converges to $T$ (up to the [[Orthogonal iteration|rotation in the complex plane]] $e^{i \Theta_k}$). So the eigenvalues of $A$ appear along the [[Schur decomposition|diagonal]] of $T_k$.

$T_k$ is initially a dense matrix. As the algorithm converges, it becomes increasingly [[Schur decomposition|upper triangular.]] The eigenvalues appear on the diagonal ordered from largest to smallest.

Note that we only converge to the [[Schur decomposition]], not the [[Diagonalizable matrices|eigendecomposition]].