# The Conjugate Gradient (CG) Method

## The Search for Coefficients and a Flawed First Attempt

Our goal is to build our solution $x_k$ from an expanding basis $\{p_1, \dots, p_k\}$ that spans the Krylov subspace ${\mathcal K}_k$. We assume $x_0 = 0$ and the Krylov process is built from $q_1 = b / \|b\|_2$.

We are looking for a basis $\{p_k\}$ such that 

$$\text{span}(p_1, \dots, p_k) = {\mathcal K}_k(A, b).$$

If we can find this basis, and also find the corresponding scalar coefficients $\mu_k$, our iterative solution $x_k$ can be built up one step at a time. The full solution $x$ would be:

$$
x = \sum_{k=1}^n \mu_k \, p_k
$$

And our $k$-th iterate is simply the partial sum:

$$
x_k = \sum_{l=1}^k \mu_l \, p_l
$$

This framework provides a beautifully simple update rule:

$$
x_{k+1} = x_k + \mu_{k+1} \, p_{k+1}
$$

The entire problem is now reduced to finding a computationally feasible method to calculate the **search directions** $p_k$ and the **step lengths** (or coefficients) $\mu_k$.

**Attempt 1: A Standard Orthogonal Basis**

Let's explore the most natural first choice for our basis: a set of **orthogonal** vectors. As we know, the Lanczos process provides exactly such a basis, $\{q_1, \dots, q_n\}$.

Let's set our search directions $p_k = q_k$.

Let $P$ be the $n \times n$ matrix whose columns are these basis vectors, $P = [p_1, \dots, p_n]$. Since we have chosen an orthogonal basis, we have $P^T P = I$.

The full solution $x$ can be written in this basis as:

$$
x = P \mu
$$

where $\mu$ is the $n \times 1$ vector of coefficients $[\mu_1, \dots, \mu_n]^T$.

Finding the coefficients $\mu$ *seems* easy. We can simply pre-multiply by $P^T$:

$$
P^T x = P^T P \mu = I \mu = \mu
$$

This gives us an exact, closed-form formula for our coefficients: $\mu_k = p_k^T x$.

But here we hit the exact same wall we encountered in the introduction. **We cannot compute $p_k^T x$ because we do not know $x$!**

This first, natural attempt has failed. It requires the very solution we are trying to find. This tells us that if our simple update $x_{k+1} = x_k + \mu_{k+1} p_{k+1}$ is to work, we need a different, more clever choice for our basis $\{p_k\}$. Simple orthogonality is not enough.

**Attempt 2: A New Optimality Condition**

Our first attempt failed because computing the coefficients $\mu = P^T x$ required knowing the solution $x$. However, we have another equation at our disposal.

We know $Ax = b$. Let's try pre-multiplying by $P^T$:

$$
P^T A x = P^T b
$$

This is the true starting point for the entire Conjugate Gradient algorithm! We know the vector $b$. Therefore, assuming we can construct the basis matrix $P$, the right-hand side $P^T b$ is a quantity we can actually compute.

Now, let's see what this means for our coefficients $\mu$. We again start with the expansion of our solution $x = P \mu$ and insert it into this new equation:

$$
P^T b = P^T (A x) = P^T A (P \mu)
$$

This gives us a $n \times n$ linear system for the unknown coefficients $\mu$:

$$
(P^T A P) \mu = P^T b
$$

We have successfully eliminated the unknown solution $x$ from the calculation. However, we are now faced with solving a new linear system. To find $\mu$, we must deal with the matrix $P^T A P$.

This brings us to our second, much more clever attempt. In Attempt 1, we chose $p_k = q_k$ (an orthogonal basis) and found that $P^T P = I$. But this didn't help us.

In Attempt 2, we will make a different choice. What if we choose our basis vectors $p_k$ such that the new system matrix, $P^T A P$, becomes trivial to solve?

The ideal scenario is to choose the basis $\{p_k\}$ such that:

$$
P^T A P = D
$$

where $D$ is a **diagonal** matrix.

This property, as we noted earlier, is called **A-orthogonality** or **conjugacy**. The basis vectors $p_k$ are "conjugate" with respect to $A$.

$$
p_k^T A p_j = 0 \quad \text{for all } k \neq j
$$

The solution for the *full* coefficient vector $\mu$ is then simply:

$$
\mu_k = \frac{(P^T b)_k}{d_k} = \frac{p_k^T b}{p_k^T A p_k}
$$

This is a breakthrough! We have found a formula for each coefficient $\mu_k$ that is entirely computable, assuming we can find the basis vectors $p_k$.

## The Meaning of $P^T A P = D$: The A-Inner Product

Let's pause to understand the condition $P^T A P = D$ more deeply. When we say that this matrix is diagonal, what are we really requiring of our basis vectors $p_i$?

This equation is a compact way of stating a new orthogonality condition. The $(i, j)$ entry of $P^T A P$ is $p_i^T A p_j$. We are therefore demanding that:

$$
p_i^T A p_j = 0 \quad \text{for all } i \neq j
$$

This leads us to define a new **inner product** (or dot product) based on our matrix $A$.

```{prf:definition} The A-Inner Product
:label: a_inner_prod

Given a Symmetric Positive Definite (SPD) matrix $A$, the **A-inner product** between two vectors $y$ and $z$ is defined as:

$$
\langle y, z \rangle_A = y^T A z
$$

This is a valid inner product precisely *because* $A$ is SPD.
```

### Geometric Interpretation of the A-Inner Product

This is not just a mathematical abstraction. It has a clear geometric meaning. Recall that since $A$ is SPD, it has an eigendecomposition:

$$
A = Q \Lambda Q^T
$$

where $Q$ is an orthogonal matrix (its columns are the eigenvectors of $A$) and $\Lambda$ is a diagonal matrix of positive eigenvalues $\lambda_i > 0$.

Let's plug this into our new inner product:

$$
\langle y, z \rangle_A = y^T A z = y^T (Q \Lambda Q^T) z = (y^T Q) \Lambda (Q^T z) = (Q^T y)^T \Lambda (Q^T z)
$$

We can group this differently using the symmetric square root of $A$, $A^{1/2} = Q \Lambda^{1/2} Q^T$:

$$
\langle y, z \rangle_A = y^T (A^{1/2} A^{1/2}) z = (A^{1/2} y)^T (A^{1/2} z)
$$

This shows that the A-inner product of $y$ and $z$ is just the **standard Euclidean dot product** of two *transformed* vectors, $A^{1/2}y$ and $A^{1/2}z$.

This transformation $v \to A^{1/2} v = (\Lambda^{1/2} Q^T) v$ can be understood in steps:

1.  **Rotate ($Q^T v$):** First, we apply $Q^T$ to the vectors, which is like rotating our coordinate system to align with the eigenvectors of $A$.
2.  **Rescale ($\Lambda^{1/2} ...$):** Second, we multiply by the diagonal matrix $\Lambda^{1/2}$, which rescales the vectors along these new axes by $\sqrt{\lambda_i}$.
3.  **Dot Product:** Finally, we take the standard dot product of these transformed vectors.

So, when we require $P^T A P = D$, we are simply saying:
**The basis vectors $p_i$ must be orthogonal with respect to this new, A-defined inner product.**

This is a very natural choice, as it's the geometry defined by the operator $A$ itself.

### The A-Norm

Just as the standard dot product defines the standard 2-norm, this new A-inner product defines a new norm, called the **A-norm**:

$$
\|z\|_A = \sqrt{\langle z, z \rangle_A} = \sqrt{z^T A z}
$$

Using our derivation from above, we can see this is equivalent to:

$$
\|z\|_A = \| \Lambda^{1/2} Q^T z \|_2 = \| A^{1/2} z \|_2
$$

As we will see, minimizing the **A-norm of the error** is the central, underlying principle of the Conjugate Gradient method.

```{note} 
We could require our basis to be fully **A-orthonormal**, meaning $p_i^T A p_j = \delta_{ij}$ (i.e., $P^T A P = I$). This would mean $\langle p_i, p_i \rangle_A = \|p_i\|_A^2 = 1$. However, for computational reasons, the standard CG algorithm uses a different normalization, which we will see shortly.
```

## Summary: The CG Algorithm Framework

Our entire derivation so far has established the core algebraic requirements for the Conjugate Gradient algorithm.

The algorithm's goal is to build a sequence of **A-orthogonal search directions** $p_i$, which are stored as the columns of a matrix $P$. This A-orthogonality is defined by the condition:

$$
P^T A P = D
$$

where $D$ is a diagonal matrix with entries $d_k = p_k^T A p_k$.

If we can find such a basis, the exact solution $x = P \mu$ can be found without knowing $x$. We compute the coefficients $\mu_k$ by solving the diagonal system:

$$
P^T b = D \mu \quad \implies \quad \mu_k = \frac{p_k^T b}{d_k}
$$

Crucially, this structure allows us to build the final solution iteratively. The $k$-th iterate $x_k$ is the optimal solution in the $k$-th subspace, and we can compute the next iterate with a simple update:

$$
x_{k+1} = x_k + \mu_{k+1} \, p_{k+1}
$$

The remaining challenge is to find an efficient way to generate this $p_k$ sequence.

## CG Optimality: A-Norm Projection

We can further interpret the solution $x_k$ in a least-squares sense.

The core principle of the CG method is that the error at step $k$, $e^{(k)} = x - x_k$, is **A-orthogonal** to the entire subspace ${\mathcal K}_k$. This is the very definition of an orthogonal projection in the $A$-inner product.

We can verify this algebraically. Using $x_k = P_k \mu^{(k)}$, we show that $P_k^T A e^{(k)} = 0$:

$$
P_k^T A (x - x_k) = P_k^T A (x - P_k \mu^{(k)}) = P_k^T b - (P_k^T A P_k) \mu^{(k)} = P_k^T b - D_k \mu^{(k)}
$$

Since $\mu^{(k)} = D_k^{-1} P_k^T b$, this becomes:

$$
P_k^T b - D_k (D_k^{-1} P_k^T b) = P_k^T b - P_k^T b = 0
$$

This result confirms that $x_k$ is the solution to a least-squares problem using the A-norm:

$$
\mu^{(k)} = \underset{y \in \mathbb{R}^k}{\rm argmin} \, \| P_k y - x \|_A, \quad \text{and} \quad x_k = P_k \mu^{(k)}
$$

This is the central optimality guarantee of the method: **CG produces the approximation in the Krylov subspace ${\mathcal K}_k$ that is closest to the true solution $x$ in the $A$-norm.**


## Computing the Search Directions $p_k$

We have established the *properties* of the A-orthogonal basis $\{p_k\}$ and *why* it gives us an optimal solution $x_k$ at every step. The final and most important piece of the puzzle is *how* to compute this sequence of vectors efficiently.

In principle, one could use a Gram-Schmidt-like process on the Krylov basis $\{b, Ab, A^2b, \ldots\}$ to make it A-orthogonal. However, this would require storing all previous vectors, which is not feasible for large systems.

A much more efficient approach exists, and it relies on the **residual vector**:

$$
r_k = b - A x_k, \quad \text{with } r_0 = b
$$

since $x_0 = 0$.

The residual is the key. Let's establish its relationship to the Krylov subspace ${\mathcal K}_k$.

Recall that ${\mathcal K}_k = \text{span}(q_1, A q_1, \dots, A^{k-1} q_1)$, where $q_1 = b/\|b\|_2$. By definition, our $k-1$ iterate $x_{k-1} \in {\mathcal K}_{k-1}$. From the properties of a Krylov subspace, if $x_{k-1} \in {\mathcal K}_{k-1}$, then $A x_{k-1} \in {\mathcal K}_k$.

Now, let's look at the residual $r_{k-1}$:

$$
r_{k-1} = b - A x_{k-1}
$$

Since $b \in {\mathcal K}_k$ (it's proportional to $q_1$) and $A x_{k-1} \in {\mathcal K}_k$, their linear combination $r_{k-1}$ must also be in ${\mathcal K}_k$.

This logic applies to all previous residuals, so we have:

$$
r_0, r_1, \dots, r_{k-1} \in {\mathcal K}_k
$$

We also know that ${\mathcal K}_k = \text{span}(p_1, \dots, p_k)$. Therefore, we have proved the following relationship:

$$
\text{span}(r_0, \dots, r_{k-1}) \subset \text{span}(p_1, \dots, p_k)
$$

This connects the residuals (which we can compute) to the search directions (which we need). We will show below that, in fact, these two subspaces are equal. This equivalence is what makes the CG algorithm so efficient.

## Key Result: Residual Orthogonality

We now prove a fundamental property of the Conjugate Gradient algorithm: the residual $r_k$ is orthogonal to the entire subspace ${\mathcal K}_k$ that was used to generate the iterate $x_k$. This is a crucial result that will lead to the algorithm's short recurrences.

```{prf:theorem} Orthogonality of Residuals
:label: cg_residual_orthogonality

The residual $r_k = b - A x_k$ is orthogonal to the subspace ${\mathcal K}_k$.

$$
r_k \perp {\mathcal K}_k
$$
```

```{prf:proof}
To prove this, we must show that $r_k$ is orthogonal to every basis vector of ${\mathcal K}_k$. Since ${\mathcal K}_k = \text{span}(p_1, \dots, p_k)$, this is equivalent to showing that $p_l^T r_k = 0$ for all $l \le k$.

Let's take the inner product for any $l$ where $1 \le l \le k$:

$$
p_l^T r_k = p_l^T (b - A x_k) = p_l^T b - p_l^T A x_k
$$

Now, we substitute the expansion of $x_k = \sum_{i=1}^k \mu_i p_i$:

$$
p_l^T r_k = p_l^T b - p_l^T A \left( \sum_{i=1}^k \mu_i p_i \right)
= p_l^T b - \sum_{i=1}^k \mu_i (p_l^T A p_i)
$$

From our framework, we know two things:
1.  $p_l^T b = \mu_l d_l$, where $d_l = p_l^T A p_l$.
2.  The basis vectors are A-orthogonal, so $p_l^T A p_i = 0$ for all $i \neq l$.

Because of the A-orthogonality, the only term that survives in the summation is when $i = l$. This gives:

$$
p_l^T r_k = (p_l^T b) - \mu_l (p_l^T A p_l)
= (d_l \mu_l) - \mu_l (d_l)
= 0
$$

Since this holds for all $l \le k$, the residual $r_k$ is orthogonal to all basis vectors of ${\mathcal K}_k$ and thus to the entire subspace.
```

## The Equivalence of Residual and Search Subspaces

We now prove a critical result that provides the "engine" for the CG algorithm's efficiency: the space spanned by the residuals is the same as the Krylov subspace.

```{prf:theorem} Residuals as a Basis
:label: cg_residual_basis

The subspace spanned by the first $k$ residuals is the same as the $k$-th Krylov subspace, which is spanned by the first $k$ search directions.

$$
\text{span}(r_0, \dots, r_{k-1}) = {\mathcal K}_k = \text{span}(p_1, \dots, p_k)
$$
```

```{prf:proof}
We already established the subset relation:

$$
\text{span}(r_0, \dots, r_{k-1}) \subset {\mathcal K}_k
$$

To prove that these two spaces are equal, we only need to show that their dimensions are equal. By definition, $\dim({\mathcal K}_k) = \dim(\text{span}(p_1, \dots, p_k)) = k$.

Thus, our goal is to prove that $\dim(\text{span}(r_0, \dots, r_{k-1})) = k$. This is true if and only if the set $\{r_0, \dots, r_{k-1}\}$ is linearly independent.

Let's use the key orthogonality result from the previous section ({prf:ref}`cg_residual_orthogonality`), which stated that $r_k \perp {\mathcal K}_k$.

* We know that $r_l \in {\mathcal K}_{l+1}$ for any $l$.
* This means that for any $l < k$, $r_l \in {\mathcal K}_{l+1} \subset {\mathcal K}_k$.
* Since $r_k \perp {\mathcal K}_k$, it must be orthogonal to all vectors within it, including all previous residuals.

Therefore, we have:

$$
r_k^T r_l = 0 \quad \text{for all } l < k
$$

This proves that the set of residual vectors $\{r_0, \dots, r_{k-1}\}$ is an **orthogonal set**.

An orthogonal set of vectors is linearly independent if and only if it does not contain the zero vector. So, we only need to show that $r_k \neq 0$.

Let's assume $r_k = 0$ for some $k \le n$.

* By definition, $r_k = b - A x_k$.
* If $r_k = 0$, then $b - A x_k = 0$, which means $A x_k = b$.
* This implies that $x_k = x$. The algorithm has converged and found the **exact solution**.

This gives us our conclusion: **Unless the algorithm has converged to the exact solution, the residual $r_k$ is non-zero.**

This confirms that $\dim(\text{span}(r_0, \dots, r_{k-1})) = k$.

Since we have a $k$-dimensional subspace $\text{span}(r_0, \dots, r_{k-1})$ that is a subset of the $k$-dimensional subspace ${\mathcal K}_k$, the two subspaces must be identical.
```

```{note}
By construction, the CG algorithm will find the exact solution in at most $n$ steps (at which point ${\mathcal K}_n = \mathbb{R}^n$). Until convergence, the set $\{r_0, \dots, r_{k-1}\}$ is an orthogonal, linearly independent set of $k$ non-zero vectors.
```

## The Three-Term Recurrence

We are now ready to derive the final, crucial component of the CG algorithm: an efficient way to compute the next search direction $p_{k+1}$. We will show that it can be computed using only the current residual $r_k$ and the previous search direction $p_k$.

We begin with two facts we have already established:

1.  **Equal Subspaces ({prf:ref}`cg_residual_basis`):** ${\mathcal K}_{k+1} = \text{span}(r_0, \dots, r_k) = \text{span}(p_1, \dots, p_{k+1})$.
2.  **Residual Orthogonality ({prf:ref}`cg_residual_orthogonality`):** $r_k \perp {\mathcal K}_k$.

From fact #1, we know that $r_k \in {\mathcal K}_{k+1}$. This means we can express $r_k$ as a linear combination of the A-orthogonal basis vectors $\{p_l\}$:

$$
r_k = \sum_{l=1}^{k+1} u_{l,k+1} \, p_l
$$

Our goal is to prove that almost all of the coefficients $u_{l,k+1}$ are zero.

Let's find a formula for these coefficients. We can isolate $u_{i,k+1}$ by taking the A-inner product of the entire equation with $p_i$ (for any $i \le k+1$):

$$
p_i^T A r_k = p_i^T A \left( \sum_{l=1}^{k+1} u_{l,k+1} \, p_l \right) = \sum_{l=1}^{k+1} u_{l,k+1} \, (p_i^T A p_l)
$$

Because the $p$-vectors are A-orthogonal, the term $p_i^T A p_l$ is zero for all $i \neq l$. The only term that survives the summation is when $l=i$, which gives $p_i^T A p_i = d_i$.

$$
p_i^T A r_k = u_{i,k+1} \, d_i
$$

This gives us an explicit formula for the coefficients:

$$
u_{i,k+1} = \frac{p_i^T A r_k}{d_i} = \frac{(A p_i)^T r_k}{d_i}
$$

Now, let's apply our second key fact, $r_k \perp {\mathcal K}_k$.

* We know that $p_i \in {\mathcal K}_i$.
* By the definition of a Krylov subspace, $A p_i \in {\mathcal K}_{i+1}$.

Consider any coefficient $u_{i,k+1}$ where **$i < k$**.

* $A p_i \in {\mathcal K}_{i+1} \subset {\mathcal K}_k$.
* Since $A p_i$ is a vector inside ${\mathcal K}_k$, and $r_k$ is orthogonal to all vectors in ${\mathcal K}_k$, their inner product must be zero:

    $$
    (A p_i)^T r_k = 0 \quad \text{for all } i < k
    $$

* This immediately implies that $u_{i,k+1} = 0$ for all $i < k$.

The expansion of $r_k$ in the $p$-basis therefore collapses. The only coefficients that are not forced to be zero are $u_{k,k+1}$ and $u_{k+1,k+1}$.

$$
r_k = u_{k,k+1} \, p_k + u_{k+1,k+1} \, p_{k+1}
$$

This is the key recurrence. By rearranging it, we can define our next search direction $p_{k+1}$ as a simple combination of the current residual $r_k$ and the previous search direction $p_k$:

$$
p_{k+1} = \frac{1}{u_{k+1,k+1}} (r_k - u_{k,k+1} \, p_k)
$$

At this point, we have not yet chosen the normalization for $p_{k+1}$. To simplify, we choose the following normalization:

$$
u_{k+1,k+1} = 1, \quad p_{k+1} = r_k - u_{k,k+1} \, p_k.
$$

**This is the key three-term recurrence relation to update $p_{k+1}$ in CG.**

With this normalization, the $p_k$ are not normalized to have unit $A$-norm. But this normalization turns out to be computationally more efficient.

## The Residuals are Orthogonal to Each Other

```{prf:theorem} Orthogonality of Residuals to Each Other
:label: cg_residuals_orthogonal

The sequence of residual vectors $\{r_0, r_1, \ldots\}$ generated by the Conjugate Gradient algorithm is an **orthogonal set**.

$$
r_k^T r_l = 0 \quad \text{for all } k \neq l
$$
```

```{prf:proof}

We will prove $r_k^T r_l = 0$ for all $k \neq l$.

Without loss of generality, let us assume $l < k$.

We will use two facts we have already established:

1.  **From {prf:ref}`cg_residual_orthogonality`:** The residual $r_k$ is orthogonal to the entire subspace ${\mathcal K}_k$.

    $$
    r_k \perp {\mathcal K}_k = \text{span}(p_1, \dots, p_k)
    $$

2.  **From {prf:ref}`cg_residual_basis`:** The residual $r_l$ is an element of the subspace ${\mathcal K}_{l+1}$.

    $$
    r_l = b - A x_l \quad \text{where } b \in {\mathcal K}_{l+1} \text{ and } x_l \in {\mathcal K}_l \implies A x_l \in {\mathcal K}_{l+1}
    $$

    Thus, $r_l \in {\mathcal K}_{l+1}$.

Since we assumed $l < k$, we get:

$$
{\mathcal K}_{l+1} \subset {\mathcal K}_k
$$

This means $r_l$ is a vector within the subspace ${\mathcal K}_k$.

From fact #1, we know that $r_k$ is orthogonal to *all* vectors in ${\mathcal K}_k$. Since $r_l$ is in ${\mathcal K}_k$, $r_k$ must be orthogonal to $r_l$.

Thus, $r_k^T r_l = 0$.
```

## Deriving the Efficient Computational Formulas

We have done the hard work of establishing the orthogonality properties of the search directions and residuals. Now, we reap the rewards.

In this section, we simplify the formulas for the coefficients $\mu_k$ and $u_{k,k+1}$. The goal is to eliminate as many matrix-vector multiplications as possible. We will discover that the parameters governing this sophisticated algorithm depend simply on the lengths (norms) of the residual vectors.

### Simplifying the Step Size $\mu_k$

Recall our formula for the optimal step size coefficient:

$$
\mu_k = \frac{p_k^T b}{d_k}
$$

We can rewrite the numerator $p_k^T b$ in a much more convenient form.
Since $r_{k-1} = b - A x_{k-1}$, we have $b = r_{k-1} + A x_{k-1}$. Substituting this into the numerator:

$$
p_k^T b = p_k^T (r_{k-1} + A x_{k-1}) = p_k^T r_{k-1} + p_k^T A x_{k-1}
$$

The second term vanishes. Why?

* $x_{k-1}$ lies in the subspace ${\mathcal K}_{k-1}$.
* By construction, the new search direction $p_k$ is **A-orthogonal** to ${\mathcal K}_{k-1}$.
* Therefore, $p_k^T A x_{k-1} = 0$.

So, $p_k^T b = p_k^T r_{k-1}$.

We can simplify further. Recall that $p_k$ is essentially the "new" part of the residual. From the 3-term recurrence relation (reversed), we know that $p_k$ differs from $r_{k-1}$ only by a multiple of the previous direction $p_{k-1}$:

$$
p_k = r_{k-1} - u_{k-1,k}p_{k-1}
$$

Using the fact that $p_{k-1} \perp r_{k-1}$ ({prf:ref}`cg_residual_orthogonality`), the inner product simplifies:

$$
p_k^T r_{k-1} = (r_{k-1} - u_{k-1,k} p_{k-1})^T r_{k-1} = \| r_{k-1} \|_2^2
$$

Substituting this into our $\mu_k$ definition gives the final, computationally efficient formula:

$$
\mu_k = \frac{\|r_{k-1}\|_2^2}{d_k}
$$

### Simplifying the Direction Update $\tau_k$

Next, we simplify the coefficient used to generate the next search direction. Recall the recurrence relation for the residual:

$$
r_k = u_{k,k+1} \, p_k + p_{k+1}
$$

We need to calculate $u_{k,k+1}$. The formula we derived was:

$$
u_{k,k+1} = \frac{p_k^T A r_k}{d_k}
$$

We can transform the term $p_k^T A$ by looking at how the residual changes.
Since $r_k = r_{k-1} - \mu_k A p_k$, we can rearrange to find $A p_k$:

$$
A p_k = \frac{1}{\mu_k} (r_{k-1} - r_k)
$$

Now, substitute this expression for $A p_k$ into the formula for $u_{k,k+1}$:

$$
u_{k,k+1} = \frac{(A p_k)^T r_k}{d_k} = \frac{1}{d_k \mu_k} (r_{k-1} - r_k)^T r_k
$$

Expanding the dot product:

$$
(r_{k-1} - r_k)^T r_k = r_{k-1}^T r_k - r_k^T r_k
$$

By the mutual orthogonality of residuals ({prf:ref}`cg_residuals_orthogonal`), $r_{k-1}^T r_k = 0$. Thus:

$$
(r_{k-1} - r_k)^T r_k = - \| r_k \|_2^2
$$

Substituting this back, and using the identity derived in step 1 ($\mu_k d_k = \|r_{k-1}\|_2^2$):

$$
u_{k,k+1} = \frac{- \| r_k \|_2^2}{\mu_k d_k} = - \frac{\| r_k \|_2^2}{\|r_{k-1}\|_2^2}
$$

This is an amazingly simple expression! It depends *only* on the ratio of the squared norms of consecutive residuals.

We define the parameter $\tau_k$ (often denoted as $\beta_k$ in standard texts) as:

$$
\tau_k = - u_{k,k+1} = \frac{\| r_k \|_2^2}{\|r_{k-1}\|_2^2}
$$

This parameter $\tau_k$ allows us to update the search direction efficiently:

$$
p_{k+1} = r_k + \tau_k p_k
$$

## The Conjugate Gradient Algorithm

After all our derivations, we arrive at a final algorithm that is both remarkably simple and computationally efficient. The various orthogonalities we proved are all implicitly enforced by the specific choice of the coefficients $\mu_k$ and $\tau_k$.

The complete CG algorithm is as follows:

Start with $x_0 = 0, r_0 = b$.
For the first iteration, set $p_1 = r_0$.

Then iterate for $k=1, 2, \dots$

$$
\begin{align*}
\mu_k &= \frac{\|r_{k-1}\|_2^2}{p_k^T A p_k} && \text{Step length} \\
x_k &= x_{k-1} + \mu_k \, p_k && \text{Update solution} \\
r_k &= r_{k-1} - \mu_k \, A p_k && \text{Update residual} \\
\tau_k &= \frac{\| r_k \|_2^2}{\|r_{k-1}\|_2^2} && \text{Direction update factor} \\
p_{k+1} &= r_k + \tau_k \, p_k && \text{Update search direction}
\end{align*}
$$

This recurrence is the computationally most efficient implementation of the CG algorithm. Each iteration requires only one sparse matrix-vector product ($A p_k$) and a few vector operations (additions and dot products). No large, dense matrices are ever formed or stored.

The CG algorithm is one of the most powerful and efficient iterative methods ever developed, but it is critical to remember that its derivation, and thus its guarantee of convergence, **applies only to Symmetric Positive Definite (SPD) matrices**.

## Key Takeaways: What to Memorize for CG

Here is a summary of the most important concepts, subspace relations, and orthogonality conditions that define the Conjugate Gradient algorithm.

### Core Definitions

* **A-Inner Product:** $\langle y, z \rangle_A = y^T A z$
* **A-Norm:** $\|z\|_A = \sqrt{z^T A z}$
* **Residual:** $r_k = b - A x_k$

### The Krylov Subspace

The entire algorithm lives in the Krylov subspace ${\mathcal K}_k$. A key insight is that this one subspace has multiple, equivalent bases that we use:

$$
{\mathcal K}_k = \text{span}(b, Ab, \dots, A^{k-1}b) = \text{span}(p_1, \dots, p_k) = \text{span}(r_0, \dots, r_{k-1})
$$

We also know that our iterate $x_k \in {\mathcal K}_k$, and if $y \in {\mathcal K}_k$, then $Ay \in {\mathcal K}_{k+1}$.

### Key Orthogonality Relations

These three properties are the "magic" that makes CG work efficiently.

* **Search Directions ($p_k$):** The search directions are **A-orthogonal**.

    $$
    p_k^T A p_l = 0 \quad (\text{for } k \neq l)
    $$

* **Residual ($r_k$):** The current residual is **orthogonal** to the entire current subspace.

    $$
    r_k \perp {\mathcal K}_k
    $$

* **Residuals ($r_k, r_l$):** The residual vectors are **mutually orthogonal**.

    $$
    r_k^T r_l = 0 \quad (\text{for } k \neq l)
    $$

### The Central Optimality Guarantee

**A-Norm Minimization:** The CG algorithm produces the approximation $x_k$ in the Krylov subspace ${\mathcal K}_k$ that is closest to the true solution $x$ in the $A$-norm.

## Orthogonality Relations in Matrix Notation

We can formalize the key properties of the CG algorithm by defining matrices whose columns are the vectors we have generated.

* $Q = [q_1, \dots, q_n]$: The orthonormal basis vectors from the Lanczos process.
* $P = [p_1, \dots, p_n]$: The A-orthogonal search direction vectors.
* $R = [r_0, \dots, r_{n-1}]$: The mutually orthogonal residual vectors.

### The Link Between Lanczos and CG Residuals

We have shown that the Krylov subspace ${\mathcal K}_k$ has two equivalent bases:
* ${\mathcal K}_k = \text{span}(q_1, \dots, q_k)$
* ${\mathcal K}_k = \text{span}(r_0, \dots, r_{k-1})$

Let's look at the "new" vector added at step $k$:
1.  **Lanczos:** The Lanczos process constructs $q_k$ to be in ${\mathcal K}_k$ and orthonormal to all previous vectors, so $q_k \perp {\mathcal K}_{k-1}$.
2.  **CG:** We proved in {prf:ref}`cg_residual_orthogonality` that $r_{k-1} \perp {\mathcal K}_{k-1}$.

Both $q_k$ and $r_{k-1}$ are non-zero vectors in ${\mathcal K}_k$ that are orthogonal to ${\mathcal K}_{k-1}$. Since this orthogonal complement (${\mathcal K}_k \ominus {\mathcal K}_{k-1}$) is a one-dimensional subspace, the two vectors must be parallel:

$$
\text{span}(q_k) = \text{span}(r_{k-1})
$$

This provides a deep connection: the CG algorithm's residual vectors are, up to scaling, the same as the orthogonal vectors generated by the Lanczos process.

### A Summary of Matrix Properties

This connection, along with our previous theorems, leads to a beautiful and concise summary of the matrix relationships.

* **Identity Matrices:**
    * $Q^T Q = I$: The Lanczos vectors are orthonormal by definition.

* **Diagonal Matrices:**
    * $R^T R = D_R$: The residuals are mutually orthogonal ({prf:ref}`cg_residuals_orthogonal`).
    * $P^T A P = D$: The search directions are A-orthogonal (the central requirement of CG).

* **Symmetric Tri-diagonal Matrices:**
    * $Q^T A Q = T$: This is the fundamental result of the Lanczos process.
    * $R^T A R = T_R$: This is the analogous (and less obvious) property for the residuals.

```{prf:proof}
$R^TAR$ is tri-diagonal because span($q_k$) = span($r_{k-1}$) and $Q^T A Q$ is tri-diagonal.

We can also prove this result directly using the orthogonality properties of the residuals. We show that $r_i^T A r_j = 0$ if $|i-j| \ge 2$.

1.  We know that $r_j \in {\mathcal K}_{j+1}$.
2.  By the property of Krylov subspaces, $A r_j \in A({\mathcal K}_{j+1}) \subset {\mathcal K}_{j+2}$.
3.  From {prf:ref}`cg_residual_orthogonality`, we know the residual $r_i$ is orthogonal to the entire subspace ${\mathcal K}_i$ (i.e., $r_i \perp {\mathcal K}_i$).
4.  Assume, without loss of generality, that $i \ge j+2$.
5.  If $i \ge j+2$, then the subspace ${\mathcal K}_{j+2}$ is a subset of ${\mathcal K}_i$. This means $A r_j$ (which is in ${\mathcal K}_{j+2}$) is also a vector in ${\mathcal K}_i$.
6.  Since $r_i \perp {\mathcal K}_i$, it must be orthogonal to every vector in it, including $A r_j$.
7.  Therefore, $r_i^T (A r_j) = 0$ for all $i \ge j+2$. By symmetry, this holds for $j \ge i+2$ as well.

```

## Connection to the Lanczos Process

The Conjugate Gradient method is not just related to the Lanczos process; it is a direct and efficient implementation of it. We can show that the CG iterate $x_k$ can be found by implicitly solving the tridiagonal system $T_k$ generated by the Lanczos iteration.

Let's derive this important connection.

1.  **Start with CG Orthogonality:** We begin with our key result from {prf:ref}`cg_residual_orthogonality`: the residual $r_k = b - A x_k$ is orthogonal to the Krylov subspace ${\mathcal K}_k$.

2.  **Use the Lanczos Basis:** The Lanczos algorithm generates an orthonormal basis $Q_k = [q_1, \dots, q_k]$ for ${\mathcal K}_k$. The orthogonality condition $r_k \perp {\mathcal K}_k$ can therefore be written in matrix form:

    $$
    Q_k^T r_k = 0 \implies Q_k^T (b - A x_k) = 0
    $$

    This rearranges to:

    $$
    Q_k^T A x_k = Q_k^T b
    $$

3.  **Express $x_k$ in the Lanczos Basis:** Since our solution $x_k$ is in ${\mathcal K}_k$, we can express it as a linear combination of the $Q_k$ basis vectors:

    $$
    x_k = Q_k y_k
    $$

    where $y_k \in \mathbb{R}^k$ is a vector of unknown coefficients.

4.  **Substitute and Solve:** We substitute this expression for $x_k$ back into our equation from step 2:

    $$
    Q_k^T A (Q_k y_k) = Q_k^T b
    $$

    Grouping the terms, we get:

    $$
    (Q_k^T A Q_k) y_k = Q_k^T b
    $$

5.  **Identify the Components:**
    * **The Matrix:** From our previous chapter on the Lanczos process, we recognize this as the definition of the $k \times k$ symmetric tridiagonal matrix $T_k$:

        $$
        T_k = Q_k^T A Q_k
        $$

    * **The Right-Hand Side:** We defined our starting vector for Lanczos as $q_1 = b / \|b\|_2$. Therefore, $Q_k^T b$ simplifies significantly. Its first component is $q_1^T b = (b^T / \|b\|_2) b = \|b\|_2$. All other components are $q_i^T b = q_i^T (\|b\|_2 q_1) = 0$ due to orthogonality.

        $$
        Q_k^T b = \begin{bmatrix} q_1^T b \\ q_2^T b \\ \vdots \\ q_k^T b \end{bmatrix} = \begin{bmatrix} \|b\|_2 \\ 0 \\ \vdots \\ 0 \end{bmatrix} = \|b\|_2 e_1
        $$

6.  **The Final System:** By substituting these components, we arrive at the final result:

    $$
    T_k y_k = \|b\|_2 e_1
    $$

This proves that the coefficients $y_k$ for the CG solution $x_k = Q_k y_k$ can be found by solving this small, symmetric tridiagonal system.

This is a profound insight. The CG algorithm is a way to solve this tridiagonal system and "build" the solution $x_k = Q_k y_k$ iteratively and implicitly, without ever forming $T_k$ or $Q_k$.

## Why Is It Called the "Conjugate Gradient" Algorithm?

The name "Conjugate Gradient" is not arbitrary; it is a literal, two-word description of the two key properties of the algorithm. Let's break down each part.

### The "Conjugate" Part

The "Conjugate" part refers to the search directions $\{p_k\}$. As we have seen, the entire algorithm is built on the requirement that these directions are A-orthogonal.

$$
p_k^T A p_l = 0 \quad \text{for } k \neq l
$$

The general mathematical term for this relationship is that the vectors $p_k$ are **conjugate** with respect to the matrix $A$.

This is the algorithm's "secret sauce." A simpler (but far less effective) method like Steepest Descent simply moves in the direction of the new gradient at each step. By enforcing A-orthogonality, the CG algorithm ensures that when we take a step in the direction $p_k$, we find the *optimal* and *final* coefficient for that basis vector. Minimizing the error in the new direction $p_{k+1}$ does not "undo" the minimization we already performed in the previous directions $p_1, \dots, p_k$.

### The "Gradient" Part

The "Gradient" part comes from the fact that the algorithm is related to gradient-based optimization methods.

To see this, we must define the function that the CG algorithm is minimizing. As we noted in our optimality proof, CG's goal is to find the $x_k$ that minimizes the **A-norm of the error**. Let's define this as our loss function $L(y)$ for an arbitrary vector $y$:

$$
L(y) = \frac{1}{2} \|x - y\|_A^2 = \frac{1}{2} (x-y)^T A (x-y)
$$

(We add the $\frac{1}{2}$ for cosmetic reasons, to cancel the '2' that appears during differentiation).

Now, let's find the **gradient** of this loss function with respect to $y$.

$$
L(y) = \frac{1}{2} (x^T A x - 2 x^T A y + y^T A y)
$$

The gradient $\nabla L(y)$ is the vector of partial derivatives:

$$
\nabla L(y) = \frac{1}{2} (-2 A x + 2 A y) = A y - A x
$$

Since $A x = b$, this simplifies to:

$$
\nabla L(y) = A y - b = -(b - A y) = -r(y)
$$

This is a profound result: **The gradient of the error function is the negative of the residual.**

This means the residual vector $r_k = b - A x_k$ that we use in our algorithm is, up to a sign, the gradient of the error function at our current iterate $x_k$. The residual $r_k$ is the direction of "steepest ascent" for the error, so $-r_k$ is the direction of "steepest descent."

### Putting It All Together

Now we can see the whole picture by looking at the update rule for the search direction:

$$
p_{k+1} = r_k + \tau_k \, p_k
$$

This equation shows that the new search direction $p_{k+1}$ is a linear combination of two vectors:

1.  **The Gradient:** The residual $r_k$, which we just proved is the gradient (up to a sign) of the error.
2.  **The Previous Direction:** The vector $p_k$.

This is why the algorithm is called **Conjugate Gradient**. It is a gradient-based method that takes the current **gradient** ($r_k$) and *corrects* it with the previous direction ($\tau_k p_k$) to create a new direction $p_{k+1}$ that is **conjugate** (A-orthogonal) to all previous directions.

## Convergence of the CG Method

We have shown that CG is an optimal and efficient algorithm. We now turn to the most practical question: how fast does it converge?

### Finite Convergence (A Theoretical Guarantee)

In a world without floating-point errors, the Conjugate Gradient algorithm is a *direct method*, not an iterative one.

```{prf:theorem} Finite Convergence of CG
:label: cg_finite_convergence

In exact arithmetic, the CG algorithm is guaranteed to find the exact solution $x$ in at most $n$ iterations, where $n$ is the dimension of $A$.
```

```{prf:proof}

The algorithm terminates when $r_k = 0$. We proved ({prf:ref}`cg_residuals_orthogonal`) that the residuals $r_0, \dots, r_{k-1}$ are mutually orthogonal. A set of $n$ non-zero, mutually orthogonal vectors forms a basis for $\mathbb{R}^n$. Therefore, the algorithm *must* find a zero residual (the exact solution) in at most $n$ steps.
```

This property is not why we use CG. For the large, sparse systems we care about, $n$ might be in the millions or billions, so $n$ iterations is impossibly expensive. Furthermore, in floating-point arithmetic, orthogonality is quickly lost, and the $n$-step guarantee fails.

The *true* goal of CG is not to run $n$ steps, but to get an excellent approximation in a number of steps $k \ll n$. This is why we treat it as an iterative method and analyze its *rate* of convergence.

### The Rate of Convergence and the Condition Number

The rate of convergence is governed by the properties of the matrix $A$, specifically its **eigenvalue distribution**. This is most easily summarized by the **condition number** $\kappa(A)$.

For an SPD matrix $A$, the condition number is the ratio of its largest to its smallest eigenvalue:

$$
\kappa(A) = \frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}
$$

A well-conditioned matrix has $\kappa(A) \approx 1$ (all eigenvalues are similar). An ill-conditioned matrix has $\kappa(A) \gg 1$ (eigenvalues are spread over many orders of magnitude).

A standard (though often pessimistic) bound on the error of the CG algorithm is given by:

$$
\|x - x_k\|_A \leq 2 \left( \frac{\sqrt{\kappa(A)} - 1}{\sqrt{\kappa(A)} + 1} \right)^k \|x - x_0\|_A
$$

### Interpreting the Error Bound

Let's analyze the term that governs the convergence: $\rho = \frac{\sqrt{\kappa(A)} - 1}{\sqrt{\kappa(A)} + 1}$. The error is reduced by this factor at each iteration $k$.

* **Well-conditioned case ($\kappa(A) \approx 1$):** If $\kappa(A)$ is close to 1, then $\sqrt{\kappa(A)} \approx 1$, and the numerator $(\sqrt{\kappa(A)} - 1)$ is very close to 0. This makes $\rho$ very small, and the error $\rho^k$ plummets to zero in just a few iterations.

* **Ill-conditioned case ($\kappa(A) \gg 1$):** If $\kappa(A)$ is very large, then $\sqrt{\kappa(A)}$ is also large. The fraction $\rho$ will be a number very close to 1. This means the error $\rho^k$ decreases *very slowly*, and we will need many iterations to achieve a good approximation.

This bound confirms our intuition: **CG converges quickly for well-conditioned matrices and slowly for ill-conditioned matrices.**

### The Real Story: Eigenvalue Clustering

This bound is often a "worst-case" scenario. The true convergence rate is more nuanced and depends on the *clustering* of eigenvalues.

* If $A$ has only $m$ distinct eigenvalues, CG will find the exact solution in at most $m$ steps.
* More practically, if $A$ has its eigenvalues "clustered" into $m$ small groups, CG will converge very quickly (as if it were solving an $m \times m$ problem).

This is the entire motivation for **preconditioning**, which is the most important topic in practical iterative methods. The goal of a preconditioner is to transform the system $Ax=b$ into a new system $\tilde{A}\tilde{x}=\tilde{b}$ where the matrix $\tilde{A}$ has the same solution but has a condition number close to 1 or has highly clustered eigenvalues.

### The Superlinear Convergence Phenomenon

In practice, the convergence of CG is often observed to accelerate as the iterations proceed. This behavior, where the error reduction per step becomes progressively larger, is known as superlinear convergence.

This acceleration is not just "fast" convergence; it is a "rate-improving" convergence. While the term "convergence rate" is a misnomer for a direct method, the phenomenon is real when CG is used iteratively. It is explained by the connection to the Lanczos process. The Lanczos algorithm (and thus CG) has a property of finding the extreme eigenvalues (and eigenvectors) first. As the CG iteration proceeds, it effectively "finds" and eliminates the error components associated with these converged "Ritz values" (eigenvalue approximations).

Once these components are removed, the algorithm behaves as if it is solving a new problem on the remaining subspace. This new, implicit problem has a smaller effective condition number ($\kappa_{eff}$). Since the convergence rate at any given stage is governed by the current $\kappa_{eff}$, the rate itself improves as $k$ increases, leading to the observed superlinear acceleration.