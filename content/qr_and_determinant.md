
# The QR Factorization and the Determinant

This is an elegant and useful theoretical result that connects QR factorization to another fundamental matrix property.

````{prf:theorem} Determinant via QR Factorization
:label: thm:det_via_qr
For any square matrix $A \in \mathbb{R}^{n \times n}$, if $A=QR$ is its QR decomposition, then:

$$
|\det(A)| = \left| \prod_{i=1}^{n} r_{ii} \right|
$$
````

````{prf:proof}
Sketch of proof:

1.  Start with the property $\det(AB) = \det(A)\det(B)$.

2.  $\det(A) = \det(QR) = \det(Q)\det(R)$.

3.  Since $Q$ is an orthogonal matrix, its determinant is always $\det(Q) = \pm 1$. Therefore, $|\det(A)| = |\det(R)|$.

4.  Since $R$ is an upper triangular matrix, its determinant is simply the product of its diagonal entries, $\det(R) = \prod_{i=1}^{n} r_{ii}$.

5.  Combining these gives the result.

````

:::{important}
**Why it's important:** This provides a very numerically stable way to compute the determinant of a matrix. The standard cofactor expansion is computationally infeasible for large matrices, and using LU factorization can suffer from overflow/underflow issues. The QR approach avoids these problems.
:::