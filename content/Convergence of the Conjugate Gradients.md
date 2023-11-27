Like the [[Lanczos process]], the convergence of [[Conjugate Gradients algorithm|CG]] is well understood.

The convergence is very fast and, to some extent, “optimal.” We can prove that:
$$
\| x - x^{(k)} \|_2 \le 
2 \, \| r^{(0)} \|_{A^{-1}} \; \Big( \frac{1 - \kappa^{-1/2}}{1 + \kappa^{-1/2}} \Big)^k
$$
where $\kappa$ is the condition number of $A$.

The number of iterations is: $k = O(\sqrt{\kappa} \; \ln \epsilon^{-1}).$

![[2022-11-28-11-17-45.png|500]]

The matrix sizes in this benchmark are equal to 512. Below are the condition numbers for each test case.

| | Condition Number
--- | ---
well-conditioned | 101
3D Laplacian | 33k
2D Laplacian | 54k