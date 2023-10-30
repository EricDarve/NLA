- Householder transformations are great for transforming an entire column to $e_1$.
- However, it is less efficient to zero out a single entry in a matrix.
- Example:

![[2022-10-11-15-08-29.png]]

Givens transformations will be very useful when dealing with Upper Hessenberg matrices (matrices that are zero below the sub-diagonal, i.e., $i > j+1$ implies that $a_{ij} = 0$) and tri-diagonal matrices.

- Left figure: upper Hessenberg
- Right figure: tri-diagonal matrix.

![[QR using Givens transformations 2023-10-18 11.41.30.excalidraw.svg]]

2 approaches:
- Householder: one big Q transform. Fastest for dense matrices.
- Givens: many small $2 \times 2$ Q transforms. Slower for dense matrices but faster for sparse matrices and for processing a single entry at a time.

Summary:
- Householder approach: one large Q
- Givens approach: many small Qs

Let's explain this on an example.
![[2022-10-11-15-10-22.png]]

Take a vector $r$ of size 2:
$$
r = \begin{pmatrix}
x \\ y
\end{pmatrix}
$$
These are the entries we want to modify. We want to zero out the 2nd entry using an orthogonal transformation.

Denote by:
$$
c = \frac{x}{\|r\|_2}, \qquad s = \frac{y}{\|r\|_2}
$$
- We have two options for the $2 \times 2$ orthogonal transformation.
	- It can be a rotation or a reflection.
	- There is no significant difference between these 2 options.

Rotation: 
$$
\begin{pmatrix}
c & s \\ -s & c
\end{pmatrix}
$$
det $= 1$

Reflection: 
$$
\begin{pmatrix}
c & s \\ s & -c
\end{pmatrix}
$$
det $= -1$

You can check that this works:
$$
\begin{pmatrix}
c & s \\ -s & c
\end{pmatrix}
\begin{pmatrix}
x \\ y
\end{pmatrix}
=
\begin{pmatrix}
c & s \\ s & -c
\end{pmatrix}
\begin{pmatrix}
x \\ y
\end{pmatrix}
=
\begin{pmatrix}
\|r\|_2 \\ 0
\end{pmatrix}
$$
Let's apply this method to our example:
![[2022-10-11-15-16-44.png]]

We get the matrix in the desired upper triangular form.

The computational cost of zeroing out a single entry is $O(n)$. So the cost for a single column is $O(n^2)$. The total computational cost for the entire matrix is $O(n^3)$.

- If matrix $A$ is a general $m \times n$ matrix, the cost is $O(mn^2)$.
- If $A$ is tri-diagonal, the cost is $O(\min(m,n))$.
- If $A$ is upper-Hessenberg, the cost is $O(mn)$.

[[QR factorization]], [[Householder transformation]], [[QR using Householder transformations]]