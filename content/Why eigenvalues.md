[[Eigenvalues|Eigenvalues]] turn a matrix multiplication into a multiplication by a scalar:
$$
A v = \lambda v
$$
[[Diagonalizable matrices|Diagonalizable matrices]] have a full eigenvector basis and can be written in the form:
$$
A = X \Lambda X^{-1}
$$
In this form, taking the power of a matrix is a simple operation:
$$
A^k = X \Lambda^k X^{-1}
$$
We can even compute the function of a matrix by using the eigendecomposition:
$$
f(A) = X f(\Lambda) X^{-1}
$$
where $f(\Lambda)$ is a diagonal matrix with entries $f(\lambda_i)$.

Eigenvalues are very useful to study time-evolving systems. Consider
$$
\frac{dx}{dt} = M x
$$
We can formally write the solution in the form:
$$
x(t) = \exp(Mt) x_0
$$
where $\exp(Mt)$ is defined as above using functions of matrices:
$$
\exp(Mt) = X \exp(\Lambda t) X^{-1}
$$
with
$$
M = X \Lambda X^{-1}
$$
and where $\exp(\Lambda t)$ is a diagonal matrix with entries $\exp(\lambda_i t)$.

This provides an exact solution to the linear system $\dot{x} = Mx$.