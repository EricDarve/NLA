- We prove that the algorithm completes and that the [[Cholesky factorization|Cholesky]] factorization exists. 
- This decomposition is also unique.

To prove existence, we use a proof by induction. Let's start with $k=1$. We have
$$
a_{11} = e_1^T A e_1 > 0
$$
This holds because $A$ is [[Symmetric Positive Definite Matrices|SPD]]. So we can write:
$$
A = [\sqrt{a_{11}} \, ] [\sqrt{a_{11}} \, ] = LL^T
$$
Now, let's consider a matrix of size $n$ and let's assume that the Cholesky factorization exists for SPD matrices of size $n-1$.

Let's perform one step of the factorization starting from $A$. We will end up with a matrix of size $n-1$ and will be able to use the induction hypothesis.

The first step of Cholesky can be written in this equivalent form:
$$
A =
\begin{pmatrix}
a_{11} & c^T \\
c & B
\end{pmatrix}
=
\begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
a_{11} & 0 \\
0 & B - \frac{1}{a_{11}} c c^T
\end{pmatrix}
\begin{pmatrix}
1 & \frac{1}{a_{11}} c^T \\
0 & I
\end{pmatrix}
$$
Note that the first column of $L$ is equal to:
$$
l_{,1} = \frac{a_{,1}}{\sqrt{a_{11}}} =
\sqrt{a_{11}}
\begin{pmatrix}
1 \\
\frac{1}{a_{11}} c
\end{pmatrix}
$$
We are using the first column of $L$ with a scaling. This is the first step in the Cholesky factorization.

We prove that $B - (1/a_{11}) c c^T$ is SPD.

1. That matrix is symmetric.
2. We now prove that $y^T (B - (1/a_{11}) c c^T) y > 0$ for any $y \neq 0.$

Check that for any row vector $z$:
$$
\begin{pmatrix}
    1 & z \\
    0 & I
\end{pmatrix}
\begin{pmatrix}
    1 & -z \\
    0 & I
\end{pmatrix}
= I
$$
Use the equation above:
$$
\begin{pmatrix}
    1 & \frac{1}{a_{11}} c^T \\
    0 & I
\end{pmatrix}
\begin{pmatrix}
    1 & -\frac{1}{a_{11}} c^T \\
    0 & I
\end{pmatrix}
= I
$$
Define
$$
X = 
\begin{pmatrix}
    1 & -\frac{1}{a_{11}} c^T \\
    0 & I
\end{pmatrix}
$$
Recall that:
$$
A =
\begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
a_{11} & 0 \\
0 & B - \frac{1}{a_{11}} c c^T
\end{pmatrix}
\begin{pmatrix}
1 & \frac{1}{a_{11}} c^T \\
0 & I
\end{pmatrix}
$$
Using the results above:
$$
X^T A X = \begin{pmatrix}
    a_{11} & 0 \\
    0 & B - \frac{1}{a_{11}} c c^T
\end{pmatrix}
$$
Take any 
$$
x = 
\begin{pmatrix}
0 \\ y
\end{pmatrix} 
\neq 0
$$
Then:
$$
x^T X^T A X x = (X x)^T A (X x) > 0
$$
because $A$ is SPD. But also:
$$
x^T X^T A X x = y^T (B - (1/a_{11}) c c^T) y
$$
So from the results above
$$
y^T (B - (1/a_{11}) c c^T) y > 0
$$
for any $y \neq 0$. So $B - (1/a_{11}) c c^T$ is SPD.

By the induction hypothesis, $B - (1/a_{11}) c c^T$ has a Cholesky factorization since its size is $n-1$:
$$
B - (1/a_{11}) c c^T
= L_1 L_1^T
$$
We get for $A$:
$$
\begin{align}
A & =
\begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
a_{11} & 0 \\
0 & B - \frac{1}{a_{11}} c c^T
\end{pmatrix}
\begin{pmatrix}
1 & \frac{1}{a_{11}} c^T \\
0 & I
\end{pmatrix} \\
\\
& = \begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
a_{11} & 0 \\
0 & L_1 L_1^T
\end{pmatrix}
\begin{pmatrix}
1 & \frac{1}{a_{11}} c^T \\
0 & I
\end{pmatrix} \\
\\
& = \begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
\sqrt{a_{11}} & 0 \\
0 & L_1
\end{pmatrix}
\begin{pmatrix}
\sqrt{a_{11}} & 0 \\
0 & L_1^T
\end{pmatrix}
\begin{pmatrix}
1 & \frac{1}{a_{11}} c^T \\
0 & I
\end{pmatrix} \\[1em]
& = L L^T
\end{align}
$$
with
$$
L = \begin{pmatrix}
1 & 0 \\
\frac{1}{a_{11}} c & I
\end{pmatrix}
\begin{pmatrix}
\sqrt{a_{11}} & 0 \\
0 & L_1
\end{pmatrix}
=
\begin{pmatrix}
\sqrt{a_{11}} & 0 \\
\frac{1}{\sqrt{a_{11}}} c & L_1
\end{pmatrix}
$$
Matrix $L$ is lower triangular as expected. This concludes the proof. $\square$

To prove uniqueness, let's assume that we have two such factorizations:
$$
L_1 L_1^T = L_2 L_2^T
$$
This implies that
$$
L_2^{-1} L_1  = L_2^T L_1^{-T}
$$
where $^{-T}$ denotes the transpose of the inverse. The left-hand side is a lower-triangular matrix while the right-hand side is an upper-triangular matrix. Therefore both matrices are diagonal. If we match the diagonal terms on the left and right-hand sides we find that:
$$
\Big( \frac{[L_1]_{ii}}{[L_2]_{ii}} \Big)^2 = 1
$$
This implies that $|[L_1]_{ii}| = |[L_2]_{ii}|.$ But both terms are positive, so they are equal. So 
$$
L_2^{-1} L_1  = I
$$
and the matrices are equal.

[[Symmetric Positive Definite Matrices]], [[Cholesky factorization]]