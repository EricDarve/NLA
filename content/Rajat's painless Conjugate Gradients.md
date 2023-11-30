
**A geometric derivation of Conjugate Gradients**

Author: **Rajat Vadiraj Dwaraknath**

Thanks to Anjan Dwaraknath for helpful discussions.

## Problem setup

We wish to solve the linear system $A x^* = b$ where $A \in \mathbb{R}^{n \times n}$ is a [[Symmetric Positive Definite Matrices|symmetric positive definite]] matrix. We use $x^*$ to denote the exact solution.

We are working in the computational model where we have access to the matrix $A$ only through [[Motivation of iterative methods for eigenvalue computation|matrix-vector products.]] That is, we have a method to compute $Av$ for any vector $v \in \mathbb{R}^n$. We also have access to the vector $b$.

Since $A$ is symmetric positive definite, it induces an [[Conjugate Gradients Version 1|inner product]] given by $\langle u, v\rangle_A = u^T A v$. Although we don't have access to the solution $x^*$, we do have access to $b = A x^*$. So, we **can** compute the expression $v^T b = v^T A x^* = \langle v, x^* \rangle_A$ for any $v$. In other words, __we can compute A-inner products $\langle x^*, v \rangle_A$ with the solution for any $v$.__

## Solution attempt

Motivated by this, we can posit a method to solve $A x^* = b$ by working in the $A$-inner product as follows:

- Compute an $A$-orthogonal basis for $\mathbb{R}^n$ which we denote $\left\{p_1, p_2, \dots, p_n\right\}$ via [[Gram-Schmidt]] in the _$A$-inner product_. Note that for this CG derivation, we don't need these basis vectors to be normalized.
- _$A$-project_ the solution $x^*$ using this basis:
$$
x^* = \sum_{i=1}^n \frac{\langle x^*, p_i \rangle_A}{\langle p_i, p_i \rangle_A} p_i
$$

This also naturally leads to an approximation scheme by truncating the sum in the projection:
$$
x_{\color{red}{k}} := \sum_{i=1}^{\color{red}{k}} \frac{\langle x^*, p_i \rangle_A}{\langle p_i, p_i \rangle_A} p_i
$$
Therefore, the sequence of approximations $x_k$ can be interpreted as the _$A$-projection_ of the solution $x^*$ onto to the sequence of increasing [[Krylov subspace|subspaces]] given by $\text{span}(p_1, \dots, p_k)$. We can iteratively update the approximations by noticing that:
$$
\begin{aligned}x_{\color{red}{k+1}} &= \sum_{i=1}^{\color{red}{k+1}} \frac{\langle x^*, p_i \rangle_A}{\langle p_i, p_i \rangle_A} p_i\\ &= \sum_{i=1}^{\color{red}{k}}\left( \frac{\langle x^*, p_i \rangle_A}{\langle p_i, p_i \rangle_A} p_i\right) + \frac{\langle x^*, p_{k+1} \rangle_A}{\langle p_{k+1}, p_{k+1} \rangle_A} p_{k+1}\\ &= x_{\color{red}k} + \frac{\langle x^*, p_{k+1} \rangle_A}{\langle p_{k+1}, p_{k+1} \rangle_A} p_{k+1}\end{aligned}
$$
Since the approximations are projections, we can also use the [[Conjugate Gradients Version 1|variational characterization]] of projection as finding the vector in the subspace that is closest to $x^*$ in the _$A$-norm_:
$$
x_k = \text{argmin}_{x \in \text{span}(p_1, \dots, p_k)} \Vert x - x^*\Vert_A
$$
Notice that there is some freedom in the choice of $\{p_i\}$ in this method. We use this freedom in a specific way to arrive at the _Conjugate Gradients_ method.

## Connecting to Conjugate Gradients

The conjugate gradients method does exactly the above procedure, but for a very specific choice of the orthogonal basis $\left\{p_1, p_2, \dots, p_n\right\}$. Specifically, it requires that these vectors span the [[Krylov subspace|Krylov sequence]] of $A$ with starting vector $b$. More precisely, **CG chooses $\{p_i\}$ such that**
$$
\text{span}(p_1, \dots, p_k) = \mathcal{K}(A, b, k) := \mathcal{K}_k \text{ for all } k
$$ 
With this choice, the [[Conjugate Gradients Version 1|variational characterization]] of the successive approximations becomes:
$$
x_k = \text{argmin}_{x \in \mathcal{K}_k} \Vert x - x^*\Vert_A
$$
which is exactly the starting definition of CG!

What remains now is to find an efficient way of computing the _$A$-orthogonal_ basis $\{p_i\}$. It turns out that choosing the successive approximation subspaces to be the Krylov sequence allows us to compute $p_i$ using a [[Three-term recurrence|short recurrence]] by connecting to the Lanczos iteration.

## The three-term recurrence for $p_i$

To compute the _$A$-orthogonal_ basis $\{p_i\}$, we can perform [[Gram-Schmidt]] in the _$A$-inner product_ on some vectors $v_1, \dots, v_n$ that span the Krylov sequence. However, Gram-Schmidt is pretty slow since to compute $p_k$ we need to _$A$-project_ out the components of $v_k$ along $p_1, \dots, p_{k-1}$ and each projection needs $O(\text{nnz}(A))$ time since we need to compute an $A$-inner product which requires multiplication by $A$. So the total time to compute the basis $\{p_1, \dots , p_k\}$ is $O(\text{nnz}(A)k^2)$. It would be nice if we only needed to project out a few components instead of all $k$ at each step. We can achieve this with a smart choice of starting vectors $v_i$.

A good choice of starting vectors $v_i$ would be one where they are already close to being $A$-orthogonal. Putting the vectors into a matrix $V_k := [v_1, \dots, v_k]$, we want the $A$-inner product of $V_k$ with its transpose to be close to a diagonal matrix. More precisely, we want $V_k^T A V_k$ to be close to a diagonal matrix. 

We have seen such a $V_k$ in the context of the [[Lanczos process|Lanczos iteration.]] Specifically, the vectors $q_1, \dots, q_k$ generated by Lanczos span the Krylov sequence and also have the property that $Q_k^T A Q_k$ is a tridiagonal matrix. This means that 
$$
q_i^T A q_j = 0 \text{ for  } |i-j| > 1 \implies q_i \perp_A q_j \text{ for  } |i-j|>1
$$
In other words, $q_k \perp_A K_{k-2}$ for all $i$.

Therefore, we choose $v_k = q_k$ for all $k$. When performing _$A$-Gram Schmidt_ on $q_k$, [[Three-term recurrence|we only need]] to project out the component of $q_k$ along $p_{k-1}$ since $q_k$ is already $A$-orthogonal to $p_{k-2}, \dots, p_1$. Therefore, we can write the Gram-Schmidt step as follows:
$$
p_k = q_k - \frac{\langle q_k, p_{k-1}\rangle_A}{\langle p_{k-1}, p_{k-1} \rangle_A} p_{k-1}
$$
This step only requires a $O(\text{nnz}(A))$ compute since we are only doing one $A$-projection. Therefore, the total time to compute $p_1, \dots, p_k$ is $O(\text{nnz}(A)k)$ which is much faster than the previous $O(\text{nnz}(A)k^2)$ time. Note that this includes the time to run the Lanczos iteration to generate $q_k$ since that also takes $O(\text{nnz}(A)k)$ time. 

Therefore, the [[Space and time costs of CG and GMRES|total time]] to compute the approximation $x_k$ is also $O(\text{nnz}(A)k)$ since we can iteratively update the approximations to obtain $x_k$ as described before, and each of these updates also only takes $O(\text{nnz}(A))$ time. 

## Bringing in the residuals

This version of CG might seem a bit different than the usual implementation since there is no mention of the residual vectors $r_k := b - Ax_k$. We can easily bring these into the picture by noticing an important fact: **the residuals $r_k$ are scaled versions of the vectors $q_{k+1}$ generated by [[Lanczos process|Lanczos.]] More precisely, $r_k$ is parallel to $q_{k+1}$ for all $0 \leq k \leq n-1$.**

Note that there is an off-by-one on the indices between $r$ and $q$ simply because in Lanczos $q_1 = b/\Vert b \Vert_2$ but the corresponding residual is $r_0 = b - A \cdot 0 = b$.

We prove this statement by showing that the ${r_k}$ form an orthogonal basis for the Krylov sequence and then use the fact that this orthogonal basis must be unique up to scaling to get the result.

First, [[All the orthogonality relations in CG|observe]] that $r_k \in K_{k+1}$. This is because $r_k = b - A x_k$ and $b \in \mathcal{K}_1$ and 
$$
x_k \in \mathcal{K}_k \implies A x_k \in \mathcal{K}_{k+1}.
$$
Now, we can rewrite the residual as 
$$
r_k = A(x^* - x_k) = A e_k
$$
where $e_k = x^* - x_k$ is the error in the approximation $x_k$. Now, since $x_k$ is the $A$-projection of $x^*$ onto the subspace $e_k$, we know by property of the projection that the error in this projection must be $A$-orthogonal to $\mathcal{K}_k$. That is, 
$$
e_k \perp_A \mathcal{K}_k \implies A e_k \perp \mathcal{K}_k \implies r_k \perp \mathcal{K}_k
$$
Combining this with $r_k \in \mathcal{K}_{k+1}$ for all $k$ gives us:
$$
r_i \perp r_j \text{ for } j<i
$$
But the orthogonality relation is symmetric, so [[All the orthogonality relations in CG|we can extend the result to:]]
$$
r_i \perp r_j \text{ for } j \neq i
$$
Therefore, ${r_i}$ is an [[Residuals and solution increments in CG|orthogonal set of vectors and]]
$$
\text{span}(r_0, \dots , r_{k-1}) = \mathcal{K}_k.
$$
This means that $r_k$ must be parallel to $q_{k+1}$ from Lanczos since we found an orthogonal basis for the Krylov sequence, and we can show by induction that this basis is unique up to scaling factors. We have therefore shown the required statement.

If we put the residuals into a matrix as follows: $R_k := [r_0, \dots, r_{k-1}]$, the above result says that $R_k^T A R_k$ is a tridiagonal matrix. 

Now, we can instead use $v_k = r_{k-1}$ for all $k$ as our starting vectors for Gram-Schmidt in the $A$-inner product when finding $p_k$. The [[Three-term recurrence|short recurrence]] for $p_k$ can now be written as:
$$
p_k = r_{k-1} - \frac{\langle r_{k-1}, p_{k-1}\rangle_A}{\langle p_{k-1}, p_{k-1} \rangle_A} p_{k-1}
$$
Additionally, we can find $r_{k+1}$ either as $b - Ax_{k+1}$, or by using the update for $x_{k+1}$ in terms of $x_k$ to get a recurrence as well:
$$
\begin{aligned}r_{k+1} &= b - Ax_{k+1}\\ &= b - A x_k  - A \frac{\langle p_{k+1}, x^* \rangle_A}{\langle p_{k+1}, p_{k+1} \rangle_A} p_{k+1}\\ &= r_k  - \frac{\langle p_{k+1}, x^* \rangle_A}{\langle p_{k+1}, p_{k+1} \rangle_A} Ap_{k+1}\end{aligned}
$$

## Putting it all together

We can now combine the [[Three-term recurrence|short recurrence]] for $p_k$ with the iterative update for $x_k$ and $r_k$ to get a more familiar version of the [[Conjugate Gradients algorithm|Conjugate Gradients]] method (we replace $k$ with $k+1$ and also write the $A$-inner products explicitly). 
$$
\begin{align}
p_{k+1} & = r_{k} - \frac{p_k^T A r_k}{p_k^T A p_k} p_{k} \\[.2em]
x_{k+1} & = x_{k} + \frac{p_{k+1}^T b}{p_{k+1}^T A p_{k+1}} p_{k+1} \\[.2em]
r_{k+1} & = r_{k} -  \frac{p_{k+1}^T b}{p_{k+1}^T A p_{k+1}} Ap_{k+1}
\end{align}
$$
with $x_0 = 0$, $r_0 = b$, $p_1 = b$.

## Summary

- We want to solve $Ax^* = b$ by using only [[Motivation of iterative methods for eigenvalue computation|matrix-vector products.]] Since we can compute inner products with $b$, we can compute $A$-inner products with $x^*$. 
- Suppose we have an [[Some orthogonality relations in CG|A-orthogonal basis]] $p_1, \dots, p_k$ for the Krylov sequence starting with $b$.
- We can obtain approximate solution $x_k$ by $A$-projecting $x^*$ onto $\text{span}(p_1, \dots, p_k) = \mathcal{K}_k$.
- To compute the basis $p_1, \dots, p_k$, do [[Gram-Schmidt]] in the $A$-inner product with [[Residuals and solution increments in CG|residuals]] $r_k$ as the starting vectors.
- This choice leads to a [[Three-term recurrence|short recurrence]] for $p_k$ resulting in an [[Conjugate Gradients algorithm|efficient algorithm]] for computing $x_k$.
