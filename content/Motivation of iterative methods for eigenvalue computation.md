A sparse matrix is a matrix that has few entries per row. Typically, for a matrix of size $n$, the number of non-zero entries per row is $O(1)$.

How can we take advantage of the many zeros present in a sparse matrix?

**Matrix-vector products with a sparse matrix can be computationally very efficient:**
$$
x \mapsto Ax
$$
Denote by $y = Ax$. We have
$$
y_i = \sum_{j=1}^n a_{ij} \, x_j
$$

- We can skip all the entries where $a_{ij} = 0$.
- Then, the  computational cost and memory are proportional to the number of non-zero entries in $A$.
- ~~LU, QR, upper Hessenberg form, QR iteration~~
- None of these methods are applicable anymore. We need a different kind of approach.
- This leads to iterative methods.
- They are less accurate and converge more slowly than methods from [[Computing eigenvalues|the previous section]], but this is the only option when $n$ becomes very large.