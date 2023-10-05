**Definition using the Leibniz formula**
$$
\det(A) = \sum_{\sigma \in S_n} \text{sgn}(\sigma) 
a_{1,\sigma(1)} \dots a_{n,\sigma(n)}
$$
$S_n$ is the set of all permutations of the set $\{ 1,2,\dots,n \}$. The signature $\text{sgn}(\sigma)$ is $+1$ if the permutation can be obtained with an even number of transpositions (i.e., exchanges of two entries); otherwise, it is $-1$.

**Properties**

- $\det(\alpha A) = \alpha^n \det(A)$
- $n$-linear function: if we fix all the columns of A except column $i$, det($A$) is a linear function of $a_i$.
- An alternating form: when two columns are identical, the determinant is 0; if you switch two columns, the determinant changes sign.
- $A$ singular if and only if det($A$) = 0.
- $\det(AB) = \det(A) \det(B)$
- $\det(A^{-1}) = \det(A)^{-1}$
- $\det(A) = \prod_{i=1}^n a_{ii}$ if $A$ is a triangular matrix.
- $\det(A) = \det(A^T)$

**Characteristic polynomial**

$$
p_A(z) = \det(zI - A)
$$
The polynomial $p_A(z)$ has degree $n$. Its $n$ complex roots are the eigenvalues of $A$.

The **multiplicity of an eigenvalue** is the number of occurrences of the corresponding root in the complete factorization of the characteristic polynomial.

**Eigenvalues**

The determinant is equal to the product of the eigenvalues:
$$
\det(A) = \prod_{i=1}^n \lambda_i
$$
In this product, eigenvalues are repeated according to their multiplicity.

**Geometric interpretation**

For any square matrix A, we consider its columns $a_i$ and the $n$-dimensional parallelepiped formed by the vectors $a_i$. The determinant of A is defined as the signed volume of this parallelepiped.

[[Eigenvalues]], [[Orthogonal matrix and projector]]