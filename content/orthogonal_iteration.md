# The Orthogonal Iteration Algorithm

The process of finding one eigenvector, building a projector, and repeating is a great theoretical concept, but it's complex and numerically difficult to implement. This process can be simplified into a single, powerful iteration: the **orthogonal iteration**. This algorithm essentially performs the "power iteration with deflation" on all eigenvectors simultaneously.

The iteration is defined by two simple steps:

1.  **Multiply by $A$:** Apply the matrix $A$ to the current set of orthonormal basis vectors $Q_k$.

    $$
    Z = A Q_k
    $$

2.  **Re-orthogonalize:** Use the $QR$ factorization to create a new orthonormal basis $Q_{k+1}$ from the resulting vectors $Z$.

    $$
    Q_{k+1} R_{k+1} = Z
    $$

    This step is what performs the implicit deflation.

  * The first column of $Q_{k+1}$ is just the normalized $A \boldsymbol{q}_1$, which is the standard power method.
  * The second column of $Q_{k+1}$ is the normalized $A \boldsymbol{q}_2$ made orthogonal to $\boldsymbol{q}_1$.
  * The $(i+1)$-th column of $Q_{k+1}$ is the normalized $A \boldsymbol{q}_i$ made orthogonal to $\text{span}\{\boldsymbol{q}_1, \dots, \boldsymbol{q}_i\}$.

This directly corresponds to our deflation idea. The column $\boldsymbol{q}_{i+1}$ from $Q_k$ is effectively converging as if it were part of a power iteration on the deflated matrix $P_{ \{\boldsymbol{q}_1, \dots, \boldsymbol{q}_i\}^\perp } A$.

## Algorithm

Here is the complete algorithm, starting from a random orthogonal matrix $Q_0$.

```python
import numpy as np

# Let A be a square matrix of size n x n
n = A.shape[0]

# Start with a random n x n orthogonal matrix
Q, _ = np.linalg.qr(np.random.randn(n, n))

# Set the number of iterations
num_iterations = 1000 

for _ in range(num_iterations):
    # 1. Multiply by A (the "power" step)
    Z = A @ Q
    
    # 2. Re-orthogonalize (the "deflation" step)
    Q, R = np.linalg.qr(Z)
    
# After iterating, Q will approximate the Schur vectors.
# The matrix T = Q.T @ A @ Q will be approximately upper triangular.
```

### Convergence

As $k \to \infty$, the matrix $Q_k$ converges to the unitary matrix $Q$ from the **Schur decomposition**, $A = Q T Q^H$.

It is important to note that $Q$ in the Schur decomposition is not unique. If $A = Q T Q^H$, we can define a new $Q' = Q D$, where $D$ is any diagonal matrix with entries $e^{i \theta_j}$ ($|\theta_j|=1$). $Q'$ is still unitary, and 

$$
A = (QD) T (QD)^H = Q (DTD^H) Q^H.
$$ 

The new matrix $T' = DTD^H$ is still upper triangular and has the same diagonal as $T$.

Because of this, we cannot rigorously say that $Q_k$ converges to a specific $Q$. Instead, we have to state the convergence more carefully:

* $Q_k$ converges to $Q$ *up to a diagonal unitary matrix*.
* More simply, the **subspaces converge**. For each $j$ from $1$ to $n$, the subspace spanned by column $j$ of $Q_k$ converges to the subspace spanned by the $j$th Schur vector.

The angle between these two subspaces goes to 0. See below for a formal definition of the angle between subspaces.

## Proof of Convergence

Let $A$ be an $n \times n$ matrix with a Schur decomposition $A = Q T Q^H$, where $Q = [\boldsymbol{v}_1 | \dots | \boldsymbol{v}_n]$ is unitary and $T$ is upper triangular. The diagonal entries of $T$ are the eigenvalues $\lambda_i = T_{ii}$.

For this proof, we make the standard assumption that the eigenvalues have distinct magnitudes:

$$
|\lambda_1| > |\lambda_2| > \dots > |\lambda_n|. 
$$

Let 

$$
Q_k = [\boldsymbol{q}_{k,1} | \dots | \boldsymbol{q}_{k,n}]
$$

be the $k$-th iterate of the algorithm. Let 

$$
\mathcal{V}_j = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_j\}
$$

be the $j$-dimensional dominant invariant subspace (spanned by the first $j$ Schur vectors). Let 

$$
\mathcal{S}_{k,j} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,j}\}
$$

be the subspace spanned by the first $j$ columns of $Q_k$.  

**Goal:** We will prove by induction that for each $j = 1, \dots, n$, the subspace $\mathcal{S}_{k,j}$ converges to the subspace $\mathcal{V}_j$ as $k \to \infty$.

### The Inductive Proof

**Inductive Hypothesis $P(j)$:** The subspace 

$$
\mathcal{S}_{k,j} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,j}\}
$$ 

converges to the invariant subspace 

$$
\mathcal{V}_j = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_j\}.
$$

### Base Case: $P(1)$

We must show that $\mathcal{S}_{k,1} = \text{span}\{\boldsymbol{q}_{k,1}\}$ converges to $\mathcal{V}_1 = \text{span}\{\boldsymbol{v}_1\}$.

Let's analyze the first column of the iteration:

1.  $\boldsymbol{z}_1 = A \boldsymbol{q}_{k,1}$
2.  From $Q_{k+1} R_{k+1} = Z$, the first column gives $\boldsymbol{q}_{k+1,1} R_{11} = \boldsymbol{z}_1$.

Since $R_{11} = \|\boldsymbol{z}_1\|_2$ (as $\boldsymbol{q}_{k+1,1}$ is a unit vector), the iteration for the first column is:

$$
\boldsymbol{q}_{k+1, 1} = \frac{A \boldsymbol{q}_{k, 1}}{\|A \boldsymbol{q}_{k, 1}\|_2}
$$

This is precisely the **standard power iteration**. Given our assumption that $|\lambda_1| > |\lambda_j|$ for all $j > 1$, the power iteration converges to the dominant eigenvector, which is the first Schur vector $\boldsymbol{v}_1$.

Thus, $\text{span}\{\boldsymbol{q}_{k,1}\} \to \text{span}\{\boldsymbol{v}_1\}$ as $k \to \infty$. The base case $P(1)$ holds.

### Inductive Step: Assume $P(i)$ holds, prove $P(i+1)$

**Assumption $P(i)$:** We assume that $\mathcal{S}_{k,i} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,i}\}$ converges to $\mathcal{V}_i = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_i\}$.

**Goal:** We must show that $P(i+1)$ holds, i.e., $\mathcal{S}_{k,i+1} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,i+1}\}$ converges to $\mathcal{V}_{i+1} = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_{i+1}\}$.

1. Let's analyze the $(i+1)$-th column of the iteration, $\boldsymbol{q}_{k+1, i+1}$. The $QR$ factorization $Q_{k+1} R_{k+1} = Z$ is equivalent to the Gram-Schmidt process. The vector $\boldsymbol{q}_{k+1, i+1}$ is computed by taking $\boldsymbol{z}_{i+1} = A \boldsymbol{q}_{k, i+1}$ and orthogonalizing it against the *preceding* vectors $\boldsymbol{q}_{k+1, 1}, \dots, \boldsymbol{q}_{k+1, i}$.

$$
\begin{aligned}
\boldsymbol{w} &= A \boldsymbol{q}_{k, i+1} - \sum_{j=1}^{i} (\boldsymbol{q}_{k+1, j}^H (A \boldsymbol{q}_{k, i+1})) \boldsymbol{q}_{k+1, j} \\
\boldsymbol{q}_{k+1, i+1} &= \frac{\boldsymbol{w}}{\|\boldsymbol{w}\|_2}
\end{aligned}
$$

The vector $\boldsymbol{w}$ is the component of $A \boldsymbol{q}_{k, i+1}$ that is orthogonal to $\mathcal{S}_{k+1, i}$.

1. Let $P^{(i)}$ be the true projection onto the orthogonal complement of $\mathcal{V}_i$, i.e., $P^{(i)} = I - \sum_{j=1}^i \boldsymbol{v}_j \boldsymbol{v}_j^H$. Let $P_k^{(i)}$ be the projection onto the orthogonal complement of $\mathcal{S}_{k,i}$.

2. From our inductive assumption $P(i)$, we know $\mathcal{S}_{k,i} \to \mathcal{V}_i$. Since the iteration $Q_{k+1} R_{k+1} = A Q_k$ is continuous, $\mathcal{S}_{k+1, i}$ also converges to $\mathcal{V}_i$. This means the projection $P_k^{(i)} \to P^{(i)}$ as $k \to \infty$.

3. Therefore, the iteration for the $(i+1)$-th column, $\boldsymbol{q}_{k+1, i+1} \propto P_{k+1}^{(i)} (A \boldsymbol{q}_{k, i+1})$, becomes asymptotically equivalent to an iteration of the form:

$$
\boldsymbol{v}^{(k+1)} \propto P^{(i)} (A \boldsymbol{v}^{(k)})
$$

where $\boldsymbol{v}^{(k)}$ represents the vector $\boldsymbol{q}_{k, i+1}$.

This is precisely the "power iteration with $PA$" that we previously discussed. We are performing a power iteration with the matrix $P^{(i)} A$ on the subspace $\mathcal{V}_i^\perp$. And we previously established that it will converge to $\boldsymbol{v}_{i+1}$ (with eigenvalue $\lambda_{i+1}$). Therefore, $\text{span}\{\boldsymbol{q}_{k, i+1}\} \to \text{span}\{\boldsymbol{v}_{i+1}\}$.

In summary, we have shown:

* $\mathcal{S}_{k,i} \to \mathcal{V}_i$ (by assumption $P(i)$)
* $\text{span}\{\boldsymbol{q}_{k, i+1}\} \to \text{span}\{\boldsymbol{v}_{i+1}\}$ (by our new finding).
* Since $\mathcal{S}_{k,i+1} = \mathcal{S}_{k,i} \oplus \text{span}\{\boldsymbol{q}_{k,i+1}\}$ and $\mathcal{V}_{i+1} = \mathcal{V}_i \oplus \text{span}\{\boldsymbol{v}_{i+1}\}$ (due to orthogonality), the convergence of the component subspaces implies the convergence of the total subspace.

Thus, $\mathcal{S}_{k,i+1} \to \mathcal{V}_{i+1}$. This proves $P(i+1)$.

By induction, the statement $P(j)$ now holds for all $j=1, \dots, n$.

This means that for each $j$, the subspace spanned by the first $j$ columns of $Q_k$ converges to the invariant subspace spanned by the first $j$ Schur vectors. This implies that the full matrix $Q_k$ converges to the Schur vector matrix $Q$ (up to a diagonal unitary matrix, as the phase of each vector $\boldsymbol{q}_{k,j}$ is not uniquely determined).