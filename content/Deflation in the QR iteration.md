We need to introduce the concept of deflation before continuing.

The goal of the [[QR iteration]] is to produce a matrix that is upper triangular. This happens in stages. Let's assume that we have applied the QR iteration to, approximately, transform matrix $A$ to a block upper triangular matrix of the form:
$$
A = \begin{pmatrix}
A_{11} & A_{12}  \\
0 & A_{22}
\end{pmatrix}
$$
We can check that the [[Eigenvalues|eigenvalues]] of $A$ are equal to those of $A_{11}$ and $A_{22}$:
$$
\lambda(A) = \lambda(A_{11}) \cup \lambda(A_{22})
$$
So if we have transformed $A$ to this form, we can get all its eigenvalues by running the [[QR iteration]] on the **smaller matrices** $A_{11}$ and $A_{22}$. This process of progressively reducing $A$ to smaller matrices is how we obtain all the eigenvalues of $A$.