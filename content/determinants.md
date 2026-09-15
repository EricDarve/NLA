# The Determinant

The **determinant** assigns a scalar $\det(A)$ to each square matrix $A$. It tells us whether $A$ is invertible and, for real matrices, how the transformation $x \mapsto Ax$ changes volume and orientation. The notation $|A|$ is also used for the determinant; it should not be confused with the absolute value $|\det(A)|$.

Throughout this section, $A$ and $B$ are $n \times n$ matrices over $\mathbb{R}$ or $\mathbb{C}$, with $n \geq 1$.

## The Leibniz Formula

We define the determinant by the **Leibniz formula**:

$$
\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma)
\prod_{i=1}^n a_{i,\sigma(i)}.
$$

Here, $S_n$ is the set of all permutations of $\{1,\ldots,n\}$. Each product selects exactly one entry from each row and each column. The **sign** $\operatorname{sgn}(\sigma)$ is $+1$ if $\sigma$ can be obtained from the identity permutation by an even number of swaps, and $-1$ if that number is odd. Although the sequence of swaps is not unique, its parity is, so the sign is well defined.

For a $1 \times 1$ matrix, $\det([a]) = a$. For a $2 \times 2$ matrix, there are two permutations, giving

$$
\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad-bc.
$$

The Leibniz formula has $n!$ terms. It is useful for proving properties, but evaluating it directly is impractical for large matrices. In the next chapter, we will use matrix factorizations to compute determinants more efficiently.

## Basic Column Properties

Write $A = [a_1,\ldots,a_n]$, where $a_j$ denotes column $j$. We first establish the rules for changing individual columns.

### Identity Matrix

$$\det(I) = 1.$$

````{prf:proof}
In the Leibniz formula for $I$, the identity permutation contributes $1$. Every other permutation selects at least one off-diagonal entry, which is zero.
````

### Linearity in Each Column

If all other columns are held fixed, the determinant is linear in the remaining column:

$$
\begin{aligned}
&\det[a_1,\ldots,\alpha u+\beta v,\ldots,a_n] \\
&\quad = \alpha\det[a_1,\ldots,u,\ldots,a_n] \\
&\qquad {}+\beta\det[a_1,\ldots,v,\ldots,a_n].
\end{aligned}
$$

This property is called **multilinearity**: the determinant is linear in each column separately.

````{prf:proof}
Each product in the Leibniz formula contains exactly one entry from the column being changed. Replacing that column by $\alpha u+\beta v$ splits each product into two terms. Summing gives the stated identity.
````

In particular, multiplying one column by $\alpha$ multiplies the determinant by $\alpha$, and a matrix with a zero column has determinant zero. Multiplying the whole matrix by $\alpha$ scales all $n$ columns, so

$$\det(\alpha A) = \alpha^n\det(A).$$

Linearity in one column does **not** mean linearity in the whole matrix. In general, $\det(A+B) \neq \det(A)+\det(B)$. For example, $\det(I_2+I_2)=4$, whereas $\det(I_2)+\det(I_2)=2$.

### Column Swaps

Exchanging two columns multiplies the determinant by $-1$.

````{prf:proof}
Let $\tau$ be the permutation that exchanges the two column indices. After the swap, the product indexed by $\sigma$ in the Leibniz formula uses the entries that were previously indexed by $\tau\circ\sigma$. These products still run through all the original products, but their signs are reversed: composing with one swap changes the sign of a permutation. Thus the new determinant is $-\det(A)$.
````

If two columns are equal, swapping them leaves the matrix unchanged while negating its determinant. Hence $\det(A)=-\det(A)$, and $\det(A)=0$.

### Adding a Multiple of Another Column

Replacing column $a_j$ by $a_j+\alpha a_i$, where $i\neq j$, does not change the determinant.

````{prf:proof}
By linearity in column $j$, the new determinant equals $\det(A)$ plus $\alpha$ times the determinant of the matrix obtained by replacing column $j$ with $a_i$. The latter matrix has two equal columns, so its determinant is zero.
````

## Key Properties

### Transpose

$$\det(A^T) = \det(A).$$

````{prf:proof}
The Leibniz formula gives

$$
\begin{aligned}
\det(A^T)
&= \sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)
   \prod_{i=1}^n a_{\sigma(i),i} \\
&= \sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)
   \prod_{j=1}^n a_{j,\sigma^{-1}(j)}.
\end{aligned}
$$

The second line just reorders the factors using $j=\sigma(i)$. A permutation and its inverse have the same sign: reversing a sequence of swaps gives the inverse using the same number of swaps. As $\sigma$ runs through $S_n$, so does $\sigma^{-1}$. The last sum is therefore the Leibniz formula for $\det(A)$.
````

Consequently, every column rule above also holds for rows: the determinant is linear in each row, a row swap reverses its sign, and adding a multiple of one row to another leaves it unchanged.

### Triangular Matrices

For an upper or lower triangular matrix,

$$\det(A) = \prod_{i=1}^n a_{ii}.$$

````{prf:proof}
Suppose $A$ is upper triangular. In a nonzero term of the Leibniz formula, the entry selected from row $n$ must be $a_{nn}$, since all other entries in that row are zero. Column $n$ has then been used, so the entry from row $n-1$ must be $a_{n-1,n-1}$. Continuing upward shows that only the identity permutation can contribute. Its term is the product of the diagonal entries. This argument also covers a zero diagonal entry, in which case every term vanishes.

If $A$ is lower triangular, $A^T$ is upper triangular with the same diagonal. Apply the result above and $\det(A)=\det(A^T)$.
````

### Multiplicativity

$$\det(AB) = \det(A)\det(B).$$

````{prf:proof}
Column $j$ of $AB$ is $\sum_{i=1}^n b_{ij}a_i$. Expanding the determinant one column at a time gives

$$
\begin{aligned}
\det(AB)
&= \sum_{i_1,\ldots,i_n=1}^n
   \left(\prod_{j=1}^n b_{i_j,j}\right) \\
&\qquad\qquad {}\cdot\det[a_{i_1},\ldots,a_{i_n}].
\end{aligned}
$$

If two indices $i_j$ coincide, the determinant in that term has repeated columns and is zero. The remaining terms correspond to permutations $\sigma$, with $i_j=\sigma(j)$. Reordering the columns of $A$ gives

$$
\det[a_{\sigma(1)},\ldots,a_{\sigma(n)}]
= \operatorname{sgn}(\sigma)\det(A).
$$

Therefore,

$$
\begin{aligned}
\det(AB)
&= \det(A)\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)
   \prod_{j=1}^n b_{\sigma(j),j} \\
&= \det(A)\det(B).
\end{aligned}
$$

The sum is the Leibniz formula for $B^T$, and $\det(B^T)=\det(B)$. No invertibility assumption was needed.
````

### Singularity

A matrix $A$ is singular (not invertible) if and only if $\det(A)=0$.

````{prf:proof}
If $A$ is singular, its columns are linearly dependent. Thus one column can be written as a linear combination of the others. Expanding the determinant in that column gives terms with repeated columns, all of which are zero. If the column itself is zero, the determinant is zero by linearity as well.

Conversely, if $A$ is invertible, multiplicativity gives

$$1 = \det(I) = \det(AA^{-1}) = \det(A)\det(A^{-1}).$$

Hence $\det(A)\neq 0$. Equivalently, a zero determinant implies that $A$ is singular.
````

### Inverse

If $A$ is invertible, then

$$\det(A^{-1}) = \frac{1}{\det(A)}.$$

````{prf:proof}
The identity $\det(A)\det(A^{-1})=1$ was established above. Since $\det(A)\neq 0$, division gives the result.
````

## The Matrix Determinant Lemma

Let $A$ be an invertible $n \times n$ matrix, and let $U,V$ be $n \times k$ matrices. The [matrix determinant lemma](https://en.wikipedia.org/wiki/Matrix_determinant_lemma#Generalization) states that

$$
\det(A+UV^T)
= \det(A)\det(I_k+V^TA^{-1}U).
$$

Here, $I_k$ is the $k \times k$ identity matrix. The update $UV^T$ has rank at most $k$, and the determinant on the right involves a $k \times k$ matrix. The identity holds even when $A+UV^T$ is singular.

````{prf:proof}
Set $X=A^{-1}U$. Factoring out $A$ gives

$$
\det(A+UV^T)=\det(A)\det(I_n+XV^T).
$$

It remains to show that $\det(I_n+XV^T)=\det(I_k+V^TX)$. Consider the $(n+k)\times(n+k)$ matrix written in four blocks:

$$
M=\begin{pmatrix}
I_n & X \\
-V^T & I_k
\end{pmatrix}.
$$

We will compute its determinant in two ways, using only additions of multiples of rows. First, add to the bottom $k$ rows the combinations of the top $n$ rows specified by $V^T$. This gives

$$
\begin{pmatrix}
I_n & X \\
0 & I_k+V^TX
\end{pmatrix}.
$$

These row operations preserve the determinant. In a nonzero term of the Leibniz formula, the bottom $k$ rows must select the last $k$ columns, leaving the first $n$ columns for the top rows. The identity block contributes $1$, so

$$\det(M)=\det(I_k+V^TX).$$

Next, start again from $M$. Subtract from the top $n$ rows the combinations of the bottom $k$ rows specified by $X$. This gives

$$
\begin{pmatrix}
I_n+XV^T & 0 \\
-V^T & I_k
\end{pmatrix}.
$$

Again, the determinant is unchanged. Now the top $n$ rows must select the first $n$ columns, and the identity block contributes $1$. Hence

$$\det(M)=\det(I_n+XV^T).$$

Equating the two expressions for $\det(M)$ and substituting $X=A^{-1}U$ proves the lemma.
````

### Rank-One Special Case

For $k=1$, write $U=u$ and $V=v$, where $u,v$ are column vectors. The matrix $I_k+V^TA^{-1}U$ then has a single entry, $1+v^TA^{-1}u$, so the lemma becomes

$$\det(A+uv^T)=(1+v^TA^{-1}u)\det(A).$$

This is the identity used in the preceding section on the [Sherman-Morrison-Woodbury formula](sherman_morrison_woodbury.md).

## Geometric Interpretation

For a real matrix $A$, $|\det(A)|$ is the volume of the parallelepiped spanned by its columns. In two dimensions this is the area of a parallelogram; in three dimensions it is the usual volume of a parallelepiped. More generally, the transformation $x\mapsto Ax$ multiplies $n$-dimensional volumes by $|\det(A)|$.

For an invertible real matrix, the sign records orientation: a positive determinant preserves orientation, while a negative determinant reverses it. If $\det(A)=0$, the columns lie in a lower-dimensional subspace and the $n$-dimensional volume is zero. The transformation then collapses at least one direction.

The volume factor need not describe how individual lengths change. A matrix can stretch one direction and compress another while leaving volume unchanged. For complex matrices, the algebraic identities above still hold, but the determinant may be complex, so this interpretation in terms of a positive or negative sign does not apply.
