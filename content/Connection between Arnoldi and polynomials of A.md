The [[Key idea of iterative methods for eigenvalue computation|vectors]] in span($Q_k$) have a special interpretation. Any vector in the span of $Q_k$ can be written as a polynomial of $A$ times $q_1.$

$$
\begin{gather}
Q_k y = \sum_{i=1}^k y_i \, q_i = \sum_{i=1}^k z_i \, A^{i-1} q_1 \\[.5em]
Q_k y = p_{k-1}(A) \, q_1
\end{gather}
$$
Polynomials of $A$ can be used to interpret the convergence of [[Arnoldi process|Arnoldi]]. Arnoldi is a process that finds $p$ such that $p(A)$ is small in some appropriate sense.

Making $p(A)$ small may start to make sense if we go back to the [[Determinant|characteristic polynomial.]] Recall the definition:
$$
p_A(z) = \det(zI - A)
$$
Then:
$$
p_A(A) = 0
$$
We provide a proof for cases where $A$ is [[Diagonalizable matrices|diagonalizable]].

Proof: Assume that $A$ is diagonalizable.
$$
\begin{gather}
A = X \Lambda X^{-1} \\[.5em]
p(A) = X \, p(\Lambda) \, X^{-1} \\[.5em]
[p(\Lambda)]_{ii} = p(\lambda_i)
\end{gather}
$$
Since
$$
p_A(\Lambda) = 0
$$
we get
$$
p_A(A) = 0
$$
$\square$

### Making $p(A)$ small

Let's explore in what sense [[Arnoldi process|Arnoldi]] makes $p(A)$ small. Here is the optimality condition statement:

Arnoldi builds a $p_k(A)$ such that $\| p_k(A) \, q_1 \|_2$ is small.

This leads to approximate eigenvalues.

Recall that we approximate the eigenvalues of $A$ using $H_k$. This is equivalent to computing the roots of the characteristic polynomial of $H_k$:
$$
p_k(z) = \det(zI - H_k)
$$

Theorem: $p_k(z) = \det(zI - H_k)$ minimizes
$$
\| p_k(A) \, q_1 \|_2
$$
among all **monic** polynomials of degree $k$:
$$
p(z) = z^k + c_{k-1} z^{k-1} + \dots + c_0
$$

**Why is this result relevant for eigenvalue approximation?**

Here is the sketch of an argument.

$$
\begin{gather}
p_k(A) \, q_1 = X p_k(\Lambda) \, X^{-1} \, q_1 \\[.5em]
p_k(\Lambda) = 
 \begin{bmatrix} \prod_{i=1}^k (\lambda_1 - \mu_i ) \\ & \prod_{i=1}^k (\lambda_2 - \mu_i ) \\ & & \ddots \\ & & & \prod_{i=1}^k (\lambda_n - \mu_i) 
\end{bmatrix}
\end{gather}
$$
So
$$
p_k(A) \text{ small } \Leftrightarrow | {\color{blue}\lambda_p} - {\color{red}\mu_i} | \text{ small.}
$$ 
But making this statement more precise is difficult for unsymmetric $A$.