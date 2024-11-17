This derivation of the Conjugate Gradient (CG) algorithm is based on Rajat's derivation. We made it self-contained so that everything you need to know is on this page. We recommend focusing your studies on this derivation. There many ways to derive CG. Of course, in the end, all these derivations are equivalent, but they will appear superficially different.

We consider $A$ symmetric positive definite and the linear system $Ax = b$. We assume that $A$ is sparse and search for an iterative method to solve the system.

As we saw previously, we can use a Lanczos process to define the Krylov subspace and the sequence of orthogonal vectors $q_1, \dots, q_k$. For this problem, we always assume that
$$
q_1 = \frac{b}{\|b\|_2}, \quad x_0 = 0
$$

Let's denote by $p_k$ a sequence of vectors such that
$$
{\rm span}(p_1, \dots, p_k) = {\rm span}(q_1, \dots, q_k)
$$
We can expand our solution $x$ in that basis:
$$
x = \sum_{i=1}^n \mu_i \, p_i
$$
If we have a method to calculate $p_i$ and $\mu_i$, then our iterative solution update is very simple:
$$
x_{k+1} = x_k + \mu_{k+1} \, p_{k+1}
$$

**Attempt 1.** Let's assume that the vectors $p_i$ are orthogonal. In fact, in that case, we have just chosen $p_i = q_i$. We then have $P^T P = I$. Since
$$
x = P \mu
$$
we have
$$
P^T x = P^T P \mu = \mu
$$
In principle this works well but we do not know $x$! So even if we can calculate the sequence $p_i$, there is no obvious way to calculate $\mu_i$.

**Attempt 2.** However there is another equation that we can use. Replace $P^Tx$ by
$$
P^T A x = P^T b
$$
This is the starting point of the entire CG algorithm! We know $b$. So if we know $P$, we can calculate $P^T b$. From there, the entire algorithm can be derived.

Recall that $x = P\mu$. If we multiply to the left by $P^T A$, we get
$$
(P^T A) x = P^T (Ax) = P^T b = (P^T A) P \mu
$$
Although we can compute $P^T b$, we now have to deal with $P^T A P$ if we want to calculate $\mu$.

In attempt 1, we had chosen $p_i = q_i$. However, other choices are possible. In attempt 2, we choose $p_i$ such that
$$
P^T A P = D
$$
where $D$ is diagonal. We will denote by $d_i$ the diagonal entries. Note that, since $A$ is SPD, we have $d_i > 0$.

What does $P^T A P$ diagonal mean? When looking at $P^T A P$, we are looking at a special dot product that uses matrix $A$. For example, the $(i,j)$ entry of $P^T A P$ is
$$
p_i^T A p_j = \langle p_i, p_j \rangle_A
$$
This dot product has the following interpretation. Recall that $A$ is SPD. So
$$
A = Q \Lambda Q^T
$$
where $Q$ is orthogonal and $\Lambda$ is diagonal with $\lambda_i > 0$. So
$$
\langle y, z \rangle_A =
y^T A z = y^T Q \Lambda Q^T z = (Q^T y)^T \Lambda (Q^T z) = (\Lambda^{1/2} Q^T y)^T (\Lambda^{1/2} Q^T z)
$$
This dot product has three steps:

1. Multiply the vectors by $Q$. That means apply a series of reflections. This is like rotating the frame of reference.
2. Multiply by the diagonal matrix $\Lambda^{1/2}$. This is a rescaling of the axes.
3. Apply the usual dot product.

So, essentially, the main thing we are doing is applying a rescaling using $\Lambda^{1/2}$.

When we say
$$
P^T A P = D
$$
we simply mean that the sequence $p_i$ is orthogonal with respect to the dot product defined by $A$. This is a very natural choice.

Note that we could require that
$$
p_i^T A p_i = 1
$$
But, for computational reasons, another normalization of $p_i$ will be used.

This new dot product allows us to define a new norm, the $A$-norm:
$$
\|z\|_A = \sqrt{z^T A z} = \| \Lambda^{1/2} Q^T z \|_2 = \| A^{1/2} z \|_2
$$

**Summary:** the CG algorithm builds a sequence of vectors $p_i$ such that
$$
P^T A P = D
$$
The exact solution $x$ is written as:
$$
x = P \mu
$$
We calculate $\mu$ using
$$
P^T b = D \mu, \quad {\rm with } \quad d_k = p_k^T A p_k
$$
or
$$
\mu_k = \frac{p_k^T b}{d_k}
$$
Then we update the solution using
$$
x_{k+1} = x_k + \mu_{k+1} \, p_{k+1}
$$

**Least-squares problem and projection.** We can further interpret the solution in a least-squares sense. From the $A$-orthogonality of $P$, we get that $x - x_k = e^{(k)}$ is $A$-orthogonal to $K_k$.

Proof:
$$
x - x_k = 
\sum_{i>k} \mu_i \, p_i
$$
Then use the fact that the $p_i$ are $A$-orthogonal.

$\square$

So in fact we are solving a least-squares problem using the $A$-norm:
$$
y_k = {\rm argmin}_{y} \| P_k \, y - x \|_A, \quad
x_k = P_k \, y_k.
$$
where $P_k$ are the first $k$ columns of $P$.
**CG produces the approximation in the Krylov subspace $K_k$  that is closest to the true solution $x$ in the $A$-norm.**

**Computing the sequence $p_k$.** In principle, using the Lanczos process, we can compute the sequence $p_k$. However, there is a more efficient approach which uses the sequence of residual vectors:
$$
r_k = b - A x_k, \quad r_0 = b.
$$
By definition, $x_k \in K_k$. From the definition of the subspace $K_k$:
$$
K_k = {\rm span}(q_1, A q_1, \dots, A^{k-1} q_1)
$$
we have that $A x_k \in K_{k+1}$. Since by definition
$$
K_k = {\rm span}(p_1, \dots, p_k)
$$
and $r_{k-1} = b - A x_{k-1},$ we have
$$
{\rm span}(r_0, \dots, r_{k-1}) = {\rm span}(p_1, \dots, p_k)
$$

Below we prove some important results involving the residuals $r_k$.

**The residual $r_k$ is orthogonal to $K_k$.** We now prove a key result: $r_k \perp K_k$.

Proof. Assume that $l \le k$. Then:
$$
p_l^T r_k = p_l^T b - p_l^T A x_k
= d_l \mu_l - p_l^T A \sum_{i=1}^k \mu_i p_i
= d_l \mu_l - d_l \mu_l = 0.
$$
$\square$

Moreover for $l > k$:
$$
p_l^T r_k = p_l^T b - p_l^T A x_k
= d_l \, \mu_l
$$
In summary, define
$$
R = \begin{bmatrix} r_0 \; r_1 \; \cdots \; r_{k-1} \; \cdots \; r_{n-1} \end{bmatrix}
$$
We can write:
$$
P^T R = L
$$
where $L$ is lower triangular, and $l_{ij} = d_i \, \mu_i$ for $i \ge j$. (Recall that column $j$ of $R$ is $r_{j-1}$.)

Note that since
$$
{\rm span}(r_0, \dots, r_{k-1}) = {\rm span}(p_1, \dots, p_k)
$$
we also have that there exists an upper triangular matrix $U$ such that
$$
R = P U.
$$
The matrix $U$ is very important and we will come back to it later.

**Three-term recurrence.** We now prove that $r_k$ is a linear combination of $p_k$ and $p_{k+1}$. From this result, we will get a short and computationally efficient recurrence formula for $p_{k+1}$.

Proof. From $R = P U$ and $P^T A P = D$, we get:
$$
(P^T A) R = (P^T A) P U = D U
$$
Since
$$
{\rm span}(p_1, \dots, p_k) = K_k
$$
we get that $Ap_k \in K_{k+1}$. We can express this result using matrix notations:
$$
AP = PH
$$
where $H$ is upper Hessenberg. Let's consider again the matrix $(P^T A) R$ from above:
$$
P^T A R = (AP)^T R = (PH)^T R
= H^T P^T R
= H^T L = W
$$
since $P^T R = L$. $H^T$ is lower Hessenberg (that is, the entries are zero for $i+2 \le j$) and $L$ is lower triangular. So $W = H^T L$ is lower Hessenberg.

Proof. Here is a more detailed proof. Consider:
$$
w_{ij} = \sum_k h_{ki} \, l_{kj}
$$
$h_{ki} = 0$ if $i+2 \le k$ and $l_{kj} = 0$ if $k \le j-1.$ So $w_{ij} = 0$ if $i + 2 \le j.$ $W$ is lower Hessenberg.

$\square$

In conclusion, $W = DU$ is both lower Hessenberg and upper triangular. So it has only two non-zero diagonals in its upper triangular part. Since $R = PU$, we have proved that (recall that column $k+1$ of $R$ is $r_k$):
$$
r_k = u_{k,k+1} \ p_k + u_{k+1,k+1} \, p_{k+1}
$$
Moreover, since $P^T A R = DU$, we have:
$$
u_{k,k+1} = \frac{p_k^T A r_k}{d_k}
$$
We now have a short recurrence relation for $p_{k+1}$:
$$
u_{k+1,k+1} \, p_{k+1} = r_k - u_{k,k+1} \, p_k
$$
So far, we have not decided on the normalization of $p_k$ yet. We now choose the simplest normalization:
$$
u_{k+1,k+1} = 1, \quad p_{k+1} = r_k - u_{k,k+1} \, p_k
$$
The $p_k$ are not normalized to have unit $A$-norm. Instead, we choose a different normalization that is computationally more efficient. With this choice **$U$ is unit upper bi-diagonal.** This means that $u_{kk} = 1$. We will show below that $u_{k,k+1} < 0$. All the other entries are 0.

**Updating the residual vectors.** We are now almost done with the complete CG algorithm. We have formulas to update $x_{k+1}$ and $p_{k+1}$. The formula to update $r_{k+1}$ can be derived from $x_{k+1}$:
$$
r_k  = b - A x_k
$$
Recall that
$$
x_{k+1} = x_k + \mu_{k+1} \, p_{k+1}
$$
Multiply by $-A$ and simplify to get:
$$
r_{k+1} = r_k - \mu_{k+1} \, A p_{k+1}
$$

**The residual $r_k$ vectors are orthogonal to each other.** There are a few more simplifications we need to make the method as computationally efficient as possible. We have seen that $R = PU$ and $P^T R = L$. We now prove that the residuals $r_k$ are orthogonal to each other.

We have:
$$
R^T R = (PU)^T R = U^T P^T R = U^T (P^T R) = U^T L
$$
Since $U^T$ and $L$ are lower triangular matrices, $R^T R$ is also lower triangular. But $R^T R$ is symmetric (and also positive definite). So $R^T R = U^T L$ is triangular and symmetric. It must be diagonal. We proved that **the residuals $r_k$ are orthogonal to each other.**

**Final simplifications.** We can now derive the final formulas for the CG algorithm. Recall that:
$$
\mu_k = \frac{p_k^T b}{d_k}
$$
But $p_k^T b = p_k^T (b - Ax_{k-1}),$ since $x_{k-1} \in K_{k-1}$ and $p_k$ is $A$-orthogonal to $K_{k-1}$. So
$$
p_k^T b = p_k^T (b - Ax_{k-1}) = p_k^T r_{k-1}
= (r_{k-1} - u_{k-1,k} \, p_{k-1})^T r_{k-1}
= r_{k-1}^T r_{k-1}
$$
since $p_{k-1}^T r_{k-1} = 0$. So we have:
$$
\mu_k = \frac{\|r_{k-1}\|_2^2}{d_k}
$$
Similarly, we can simplify
$$
u_{k,k+1} = \frac{p_k^T A r_k}{d_k}
$$
From the previous relations, we have
$$
A p_k = \mu_k^{-1} (r_{k-1} - r_k)
$$
And $\mu_k d_k = \|r_{k-1}\|_2^2$. So using the orthogonality of $r_k$
$$
u_{k,k+1} = \frac{p_k^T A r_k}{d_k}
= \frac{(Ap_k)^T r_k}{d_k}
= \frac{(r_{k-1} - r_k)^T r_k}{\mu_k d_k}
= - \frac{\| r_k \|_2^2}{\|r_{k-1}\|_2^2} < 0
$$
This is an amazingly simple expression! Below we will denote by $\tau_k = - u_{k,k+1}$.

**The Conjugate Gradient Algorithm.** The complete CG algorithm is as follows. Start with
$$
x_0 = 0, \quad r_0 = b, \quad p_1 = b.
$$
Then iterate starting, from $k=1$:
$$
\begin{align*}
\mu_k &= \frac{\|r_{k-1}\|_2^2}{p_k^T A p_k} \\
x_k &= x_{k-1} + \mu_k \, p_k \\
r_k &= r_{k-1} - \mu_k \, A p_k \\
\tau_k &= \frac{\| r_k \|_2^2}{\|r_{k-1}\|_2^2} \\
p_{k+1} &= r_k + \tau_k \, p_k
\end{align*}
$$
This recurrence is the computationally most efficient implementation of the CG algorithm. It relies on sparse matrix vector products with $A$ and just very few vector operations. The CG algorithm is one of the most efficient iterative methods for solving linear systems. But note that it only applies to SPD matrices.

**Summary of key equations.** We list all the key results we have used in our derivation of the CG algorithm:
$$
\begin{align*}
x &= P \mu \\
P^T A P &= D, \quad \text{where $D$ is diagonal} \\
\|z\|_A &= \sqrt{z^T A z} = \| A^{1/2} z \|_2 \\
r_k &= b - A x_k \\
R &= P U, \quad \text{where $U$ is unit upper bi-diagonal} \\
P^T R &= L, \quad \text{where $L$ is lower triangular}
\end{align*}
$$
We have:
$$
K_k = {\rm span}(b, Ab, \dots, A^{k-1}b) = {\rm span}(q_1, \dots, q_k) = {\rm span}(p_1, \dots, p_k) = {\rm span}(r_0, \dots, r_{k-1}) 
$$
and $x_k \in K_k$. For any vector $y \in K_k$, $Ay \in K_{k+1}$.

Key orthogonality relations:

- The vectors $p_k$ are $A$-orthogonal.
- The residual $r_k$ is orthogonal to $K_k$.
- The residuals $r_k$ are orthogonal to each other.

Key optimality relation:

- The CG algorithm produces the approximation $x_k$ in the Krylov subspace $K_k$ that is closest to the true solution $x$ in the $A$-norm.

**Why is it called the Conjugate Gradient algorithm?** The name comes from the fact that the directions $p_k$ are $A$-orthogonal or **conjugate** to each other. 

Moreover consider the loss function
$$
L(y) = \| x - y \|_A^2 = (x-y)^T A (x-y)
$$
Let's calculate the gradient with respect to $y$:
$$
\nabla L(y) = 2 A (y - x) = -2 (b - A y) = -2 r(y)
$$
Our residual $r_k = b - A x_k$ is parallel to the gradient of the loss function at $x_k$. Recall our update equation
$$
p_{k+1} = r_k + \tau_k \, p_k
$$
The new direction $p_{k+1}$ is the optimal direction to update $x_{k+1}$. This equation means that $p_{k+1}$ is a linear combination of the gradient at $x_k$ and the previous direction $p_k.$

This is why the CG algorithm is called the Conjugate Gradient algorithm.