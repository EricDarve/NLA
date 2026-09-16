# Block Matrix Operations

A **block matrix** is an ordinary matrix whose rows and columns have been grouped into smaller matrices, called **blocks**. The partition changes how we write the matrix, but not its entries or its action on a vector.

For example,

$$
A=\left(\begin{array}{cc|c}
1&2&3\\
4&5&6\\ \hline
7&8&10
\end{array}\right)
=\begin{pmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{pmatrix}.
$$

Here $A_{11}$ is $2\times2$, $A_{12}$ is $2\times1$, $A_{21}$ is $1\times2$, and $A_{22}$ is $1\times1$. Blocks need not be square or have the same size. However, blocks in the same block row must have the same number of rows, and blocks in the same block column must have the same number of columns.

More generally, if the row groups have sizes $m_1,\ldots,m_p$ and the column groups have sizes $n_1,\ldots,n_q$, then

$$A_{ij}\in\mathbb{R}^{m_i\times n_j}.$$

All the rules below also apply to complex matrices.

## Addition, Scalar Multiplication, and Transpose

If $A$ and $B$ have the same size and the same partition, addition and scalar multiplication act on each block:

$$
\begin{aligned}
(A+B)_{ij}&=A_{ij}+B_{ij},\\
(cA)_{ij}&=cA_{ij}.
\end{aligned}
$$

Transposing a matrix exchanges its block rows and block columns **and transposes each block**:

$$
A^T=\begin{pmatrix}
A_{11}^T&A_{21}^T\\
A_{12}^T&A_{22}^T
\end{pmatrix}.
$$

The same rule holds for the conjugate transpose, with $T$ replaced by $H$.

## Block Matrix Multiplication

For $AB$, the column groups of $A$ must match the row groups of $B$. Specifically, if

$$
A_{ik}\in\mathbb{R}^{m_i\times n_k},\qquad
B_{kj}\in\mathbb{R}^{n_k\times \ell_j},
$$

then each product $A_{ik}B_{kj}$ has size $m_i\times\ell_j$, and

$$ (AB)_{ij}=\sum_{k=1}^q A_{ik}B_{kj}. $$

This is ordinary matrix multiplication with its sums grouped according to the partition. The factors must stay in their original order: blocks are matrices, so they generally do not commute.

For two matrices with compatible $2\times2$ block partitions, write $C=AB$. Its four blocks are

$$
\begin{aligned}
C_{11}&=A_{11}B_{11}+A_{12}B_{21},\\
C_{12}&=A_{11}B_{12}+A_{12}B_{22},\\
C_{21}&=A_{21}B_{11}+A_{22}B_{21},\\
C_{22}&=A_{21}B_{12}+A_{22}B_{22}.
\end{aligned}
$$

### Multiplying a Block Vector

Partitioning a vector to match the columns of $A$ gives

$$
A\begin{pmatrix}x_1\\x_2\end{pmatrix}
=\begin{pmatrix}
A_{11}x_1+A_{12}x_2\\
A_{21}x_1+A_{22}x_2
\end{pmatrix}.
$$

For the matrix at the start of this section, take $x_1=(1,-1)^T$ and $x_2=2$. Then

$$
\begin{aligned}
A_{11}x_1+A_{12}x_2
&=\begin{pmatrix}-1\\-1\end{pmatrix}
 +\begin{pmatrix}6\\12\end{pmatrix}
 =\begin{pmatrix}5\\11\end{pmatrix},\\
A_{21}x_1+A_{22}x_2&=-1+20=19.
\end{aligned}
$$

Thus $Ax=(5,11,19)^T$, exactly as in the usual calculation.

Another useful partition groups $A$ by columns and $B$ by rows. If $a_k$ is column $k$ of $A$ and $b_k^T$ is row $k$ of $B$, then

$$AB=\sum_{k=1}^n a_kb_k^T.$$

Each outer product $a_kb_k^T$ contributes one term to the matrix product.

## Block Diagonal and Block Triangular Matrices

A **block diagonal matrix** has square diagonal blocks and zero blocks elsewhere:

$$
D=\operatorname{diag}(D_1,\ldots,D_r)
=\begin{pmatrix}
D_1&0&\cdots&0\\
0&D_2&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
0&0&\cdots&D_r
\end{pmatrix}.
$$

The blocks can have different sizes. The system $Dx=b$ separates into independent systems $D_ix_i=b_i$. Consequently, $D$ is invertible exactly when every $D_i$ is invertible, and in that case

$$D^{-1}=\operatorname{diag}(D_1^{-1},\ldots,D_r^{-1}).$$

Multiplying the two block diagonal matrices verifies this formula. Also,

$$
\det(D)=\prod_{i=1}^r\det(D_i),\qquad
\operatorname{tr}(D)=\sum_{i=1}^r\operatorname{tr}(D_i).
$$

The notation $D_1\oplus\cdots\oplus D_r$, called a **direct sum**, is another way to write this block diagonal matrix.

A **block upper triangular matrix** has zero blocks below its block diagonal. For two square diagonal blocks,

$$
T=\begin{pmatrix}A&B\\0&D\end{pmatrix},\qquad
\det(T)=\det(A)\det(D).
$$

````{prf:proof}
In the determinant expansion, a nonzero term must select the entries from the bottom rows in the columns belonging to $D$, because the bottom-left block is zero. These selections use all those columns, so the top rows must select entries from $A$. The sum therefore factors into $\det(A)\det(D)$.
````

The result extends to any number of diagonal blocks, and to block lower triangular matrices by transposing. It includes the block diagonal determinant formula above. The diagonal blocks themselves need not be triangular.

If $A$ and $D$ are invertible, the system $Tx=b$ can be solved one block at a time:

$$
\begin{aligned}
Dx_2&=b_2,\\
Ax_1&=b_1-Bx_2.
\end{aligned}
$$

We will develop this procedure further when studying triangular systems.

## The Schur Complement and Block Inversion

Consider a square matrix

$$
M=\begin{pmatrix}A&B\\C&D\end{pmatrix},
$$

where $A$ is $p\times p$, $D$ is $q\times q$, and the other blocks have the corresponding sizes. Assume $A$ is invertible. The system

$$
M\begin{pmatrix}x\\y\end{pmatrix}
=\begin{pmatrix}f\\g\end{pmatrix}
$$

consists of two equations:

$$
\begin{aligned}
Ax+By&=f,\\
Cx+Dy&=g.
\end{aligned}
$$

Solve the first for $x$ and substitute into the second:

$$
\begin{aligned}
x&=A^{-1}(f-By),\\
(D-CA^{-1}B)y&=g-CA^{-1}f.
\end{aligned}
$$

The matrix

$$S=D-CA^{-1}B$$

is the **Schur complement of $A$ in $M$**. It describes the remaining system after eliminating $x$.

The same elimination can be written as a block factorization:

$$
M=
\begin{pmatrix}I_p&0\\CA^{-1}&I_q\end{pmatrix}
\begin{pmatrix}A&B\\0&S\end{pmatrix}.
$$

Block multiplication verifies the factorization: the bottom-right block is $CA^{-1}B+S=D$. Taking determinants gives

$$\det(M)=\det(A)\det(S).$$

Thus, under the assumption that $A$ is invertible, $M$ is invertible exactly when $S$ is invertible.

If $S$ is invertible, the inverse has the block form

$$
M^{-1}=\begin{pmatrix}E&F\\G&H\end{pmatrix},
$$

where the inverse Schur complement determines all four blocks:

$$
\begin{aligned}
H&=S^{-1}=(D-CA^{-1}B)^{-1},\\
E&=A^{-1}(A+BHC)A^{-1},\\
F&=-A^{-1}BH,\\
G&=-HCA^{-1}.
\end{aligned}
$$

Here $H$ inverts the reduced system, while $E$ is $A^{-1}$ plus a correction from the coupling blocks $B$ and $C$. The off-diagonal formulas have matching forms. In particular, if $M$ is real symmetric, then $E$ and $H$ are symmetric and $G=F^T$. The block inverse can be verified directly by multiplication.

These formulas require both $A$ and $S$ to be invertible. An invertible matrix need not have invertible diagonal blocks: for example, $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ is its own inverse. Also, inverting a general block matrix does not mean inverting its blocks separately.

In computations, expressions such as $A^{-1}B$ are evaluated by solving $AX=B$. The role of block elimination in [LU factorization](lu_decomposition.md) will be explained later in the course.
