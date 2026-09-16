# The Four Fundamental Subspaces

Let $A\in\mathbb{R}^{m\times n}$ represent a linear map from $\mathbb{R}^n$ to $\mathbb{R}^m$. Its four fundamental subspaces describe which outputs are possible, which inputs map to zero, and how these spaces fit together. Two lie in the input space $\mathbb{R}^n$, and two lie in the output space $\mathbb{R}^m$.

We use $R(A)$ for the range and $N(A)$ for the null space. These are also written $\operatorname{range}(A)$ and $\ker(A)$.

## 1. The Column Space (Range)

The **column space**, or **range**, consists of all possible outputs:

$$R(A)=\{Ax:x\in\mathbb{R}^n\}\subseteq\mathbb{R}^m.$$

If the columns of $A$ are $a_1,\ldots,a_n$, then

$$Ax=x_1a_1+\cdots+x_na_n.$$

Thus $R(A)=\operatorname{span}(a_1,\ldots,a_n)$. Its dimension is the **rank** $r$ of $A$: the maximum number of linearly independent columns.

## 2. The Null Space (Kernel)

The **null space**, or **kernel**, consists of all inputs mapped to zero:

$$N(A)=\{x\in\mathbb{R}^n:Ax=0\}.$$

Its dimension is called the **nullity**. A vector in $N(A)$ can be added to any input without changing its output. More precisely,

$$Ax=Ay\quad\Longleftrightarrow\quad x-y\in N(A).$$

Therefore, the map is one-to-one exactly when $N(A)=\{0\}$.

## 3. The Row Space

The **row space** is the span of the rows of $A$, written as column vectors. Equivalently, it is the range of the transpose:

$$R(A^T)=\{A^Ty:y\in\mathbb{R}^m\}\subseteq\mathbb{R}^n.$$

Its dimension is also $r$. We will prove below that row rank and column rank are equal.

## 4. The Left Null Space

The **left null space** is the null space of the transpose:

$$N(A^T)=\{y\in\mathbb{R}^m:A^Ty=0\}.$$

The name comes from the equivalent equation $y^TA=0^T$. Such a vector gives a linear combination of the rows of $A$ that equals zero.

## Orthogonal Complements

The four subspaces form two pairs of orthogonal complements:

$$
\begin{aligned}
N(A)&=R(A^T)^\perp,\\
N(A^T)&=R(A)^\perp.
\end{aligned}
$$

````{prf:proof}
The equation $Ax=0$ says that the dot product of $x$ with every row of $A$ is zero. This holds exactly when $x$ is perpendicular to every vector in the row space. Hence $N(A)=R(A^T)^\perp$.

Applying the same argument to $A^T$ gives $N(A^T)=R(A)^\perp$.
````

Being orthogonal complements means both that the two spaces are perpendicular and that together they span the entire input or output space:

$$
\begin{aligned}
\mathbb{R}^n&=R(A^T)\oplus N(A),\\
\mathbb{R}^m&=R(A)\oplus N(A^T).
\end{aligned}
$$

Here $\oplus$ denotes a direct sum. Every vector has a unique decomposition into components in the two indicated spaces, and these components are orthogonal. For example, every input can be written uniquely as

$$x=x_{\mathrm{row}}+x_{\mathrm{null}},$$

with $x_{\mathrm{row}}\in R(A^T)$ and $x_{\mathrm{null}}\in N(A)$. Then

$$Ax=Ax_{\mathrm{row}}.$$

The null-space component has no effect on the output.

## Dimensions and the Rank-Nullity Theorem

The **rank-nullity theorem** states that

$$\dim R(A)+\dim N(A)=n.$$

````{prf:proof}
Choose a basis $z_1,\ldots,z_k$ of $N(A)$ and extend it to a basis

$$z_1,\ldots,z_k,v_1,\ldots,v_{n-k}$$

of $\mathbb{R}^n$. Since $Az_i=0$, the vectors $Av_1,\ldots,Av_{n-k}$ span $R(A)$.

They are also linearly independent. Indeed, if

$$\sum_{j=1}^{n-k}c_jAv_j=0,$$

then $\sum_j c_jv_j$ belongs to $N(A)$ and is therefore a linear combination of the $z_i$. Independence of the extended basis forces every $c_j$ to be zero. Thus $\dim R(A)=n-k$, proving the result.
````

Since $R(A^T)$ and $N(A)$ are orthogonal complements, their dimensions also add to $n$. Consequently,

$$\dim R(A^T)=n-\dim N(A)=r.$$

This proves equality of row and column rank. The output-space decomposition then gives $\dim N(A^T)=m-r$.

| Subspace | Ambient space | Dimension |
|---|---|---|
| Column space $R(A)$ | $\mathbb{R}^m$ | $r$ |
| Null space $N(A)$ | $\mathbb{R}^n$ | $n-r$ |
| Row space $R(A^T)$ | $\mathbb{R}^n$ | $r$ |
| Left null space $N(A^T)$ | $\mathbb{R}^m$ | $m-r$ |

In particular, $0\le r\le\min(m,n)$. These dimension and orthogonality statements are often called the **fundamental theorem of linear algebra**.

The restriction of $A$ to its row space is a one-to-one map onto its column space. It reaches every output because $Ax=Ax_{\mathrm{row}}$, and it is one-to-one because $R(A^T)\cap N(A)=\{0\}$. This does not mean that $A$ preserves lengths or angles.

## A Rectangular Example

Consider

$$
A=\begin{pmatrix}
1&0&1&0\\
0&1&0&1\\
1&1&1&1
\end{pmatrix}.
$$

The first two columns are independent, and the last two repeat them, so $r=2$. Bases for the column and row spaces are

$$
\begin{aligned}
R(A)&=\operatorname{span}\left\{
\begin{pmatrix}1\\0\\1\end{pmatrix},
\begin{pmatrix}0\\1\\1\end{pmatrix}\right\},\\
R(A^T)&=\operatorname{span}\left\{
\begin{pmatrix}1\\0\\1\\0\end{pmatrix},
\begin{pmatrix}0\\1\\0\\1\end{pmatrix}\right\}.
\end{aligned}
$$

To find the null space, $Ax=0$ gives $x_1+x_3=0$ and $x_2+x_4=0$; the third equation is their sum. Therefore,

$$
N(A)=\operatorname{span}\left\{
\begin{pmatrix}-1\\0\\1\\0\end{pmatrix},
\begin{pmatrix}0\\-1\\0\\1\end{pmatrix}\right\}.
$$

Similarly, $A^Ty=0$ gives $y_1+y_3=y_2+y_3=0$, so

$$
N(A^T)=\operatorname{span}\left\{
\begin{pmatrix}-1\\-1\\1\end{pmatrix}\right\}.
$$

The dimensions are $2,2,2,1$, respectively. Taking dot products verifies that the displayed row-space basis is perpendicular to the null-space basis, and the column-space basis is perpendicular to the left-null-space basis.

## What the Subspaces Tell Us about $Ax=b$

A solution exists exactly when $b\in R(A)$. Equivalently,

$$y^Tb=0\quad\text{for every }y\in N(A^T).$$

This follows from $R(A)=N(A^T)^\perp$. In the example above, the condition is $b_3=b_1+b_2$.

If $x_0$ is one solution, then all solutions have the form

$$x=x_0+z,\qquad z\in N(A).$$

Thus a consistent system has a unique solution exactly when $r=n$ (**full column rank**). Every right-hand side has a solution exactly when $r=m$ (**full row rank**). For a square matrix, both conditions are equivalent to invertibility.

For each consistent right-hand side, there is exactly one solution $x_{\mathrm{row}}$ in the row space. It also has the smallest Euclidean norm among all solutions, since orthogonality gives

$$\|x_{\mathrm{row}}+z\|_2^2
=\|x_{\mathrm{row}}\|_2^2+\|z\|_2^2,
\qquad z\in N(A).$$

When $b\notin R(A)$, its orthogonal projection onto $R(A)$ is the closest attainable output. We will use this observation in the chapter on [least squares](least_squares.md).

## Complex Matrices

For $A\in\mathbb{C}^{m\times n}$, use the conjugate transpose $A^H$ and the inner product $x^Hy$. The orthogonal decompositions become

$$
\begin{aligned}
\mathbb{C}^n&=R(A^H)\oplus N(A),\\
\mathbb{C}^m&=R(A)\oplus N(A^H).
\end{aligned}
$$

Here $R(A^H)$ is spanned by the conjugate transposes of the rows of $A$, and $N(A^H)$ consists of vectors satisfying $y^HA=0$. The dimension formulas remain the same, with dimensions taken over $\mathbb{C}$.
