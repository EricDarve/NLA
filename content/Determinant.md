The determinant is equal to the product of the eigenvalues:
$$
\det(A) = \prod_{i=1}^n \lambda_i
$$
Geometric interpretation. For any square matrix A, we consider its columns $a_i$ and the $n$-dimensional parallelepiped formed by the vectors $a_i$. The determinant of A is defined as the signed volume of this parallelepiped.

Properties:
- $\det(AB) = \det(A) \det(B)$
- $\det(\alpha A) = \alpha^n \det(A)$
- An alternating form: when two columns are identical, the determinant is 0; if you switch two columns, the determinant changes sign
- $n$-linear function: if we fix all the columns of A except column $i$, det($A$) is a linear function of $a_i$.
- $A$ singular if and only if det($A$) = 0.
- $\det(A) = \det(A^T)$
- $\det(A^{-1}) = \det(A)^{-1}$
- $\det(A) = \prod_{i=1}^n a_{ii}$ if $A$ is a triangular matrix.

**Characteristic polynomial**
$$
p_A(z) = \det(zI - A)
$$
The polynomial $p_A(z)$ has degree $n$. Its $n$ complex roots are the eigenvalues of $A$.

[[Eigenvalues]], [[Orthogonal matrix and projector]]