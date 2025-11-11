# Splitting Methods and Convergence Theory

## Definition of a Splitting Method

A **stationary iterative method** for solving a linear system $A\mathbf{x} = \mathbf{b}$ begins by expressing the coefficient matrix $A$ as the difference of two matrices, $M$ and $N$. This decomposition is known as a **splitting** of $A$:

$$A = M - N$$

For the splitting to be useful for iteration, the matrix $M$ must be chosen such that it is nonsingular. $M$ is frequently referred to as the **preconditioning matrix** or sometimes $C$ in the context of splitting $A=C-R$.

This decomposition allows the original equation $A\mathbf{x} = \mathbf{b}$ to be rewritten as $M\mathbf{x} = N\mathbf{x} + \mathbf{b}$. This linear system directly defines the iterative procedure:

$$M \mathbf{x}^{(k+1)} = N \mathbf{x}^{(k)} + \mathbf{b}$$

If we multiply by $M^{-1}$ (which exists since $M$ is nonsingular), the iterative step can be explicitly written in the fixed-point form:

$$\mathbf{x}^{(k+1)} = M^{-1} N \mathbf{x}^{(k)} + M^{-1} \mathbf{b} = R \mathbf{x}^{(k)} + \mathbf{c}$$

Here, $R = M^{-1}N$ is the **iteration matrix** and $\mathbf{c} = M^{-1}\mathbf{b}$. Common iterative methods, such as Jacobi, Gauss-Seidel, and SOR, are all derived from specific choices of this splitting.

## Necessary and Sufficient Condition for Convergence

The fundamental question of whether a splitting method converges hinges entirely on the properties of the iteration matrix $R$.

```{prf:theorem} Convergence of Splitting Methods
:label: thm:splitting_convergence

These iterative schemes converge for any initial guess $x^{(0)}$ and any right-hand side vector $b$ if and only if the **spectral radius** of the iteration matrix $R = M^{-1}N$ is strictly less than one:

$$\rho(R) < 1$$

The spectral radius $\rho(R)$ is defined as the maximum absolute value of all eigenvalues $\lambda$ of $R$.
```

```{prf:proof}

The convergence proof relies on analyzing the behavior of the error vector $\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}$, where $\mathbf{x}$ is the exact solution to $A\mathbf{x}=\mathbf{b}$.

## 1. Establishing the Error Recursion

The exact solution $\mathbf{x}$ satisfies the steady-state equation derived from the splitting:

$$M \mathbf{x} = N \mathbf{x} + \mathbf{b}$$

The iterative solution satisfies the recurrence relation:

$$M \mathbf{x}^{(k+1)} = N \mathbf{x}^{(k)} + \mathbf{b}$$

Subtracting the exact equation from the iteration equation yields the relationship between consecutive error vectors:

$$
\begin{aligned}
M (\mathbf{x}^{(k+1)} - \mathbf{x}) &= N (\mathbf{x}^{(k)} - \mathbf{x})\\
M \mathbf{e}^{(k+1)} &= N \mathbf{e}^{(k)}
\end{aligned}
$$

Since $M$ is nonsingular, we can isolate $\mathbf{e}^{(k+1)}$:

$$\mathbf{e}^{(k+1)} = M^{-1} N \mathbf{e}^{(k)} = R \mathbf{e}^{(k)}$$

By repeated substitution, the error after $k$ steps depends linearly on the initial error $\mathbf{e}^{(0)}$ and the $k$-th power of the iteration matrix:

$$\mathbf{e}^{(k)} = R^k \mathbf{e}^{(0)}$$

Convergence requires that $\lim_{k\to\infty} \mathbf{e}^{(k)} = \mathbf{0}$, which necessitates that the sequence of matrix powers $R^k$ tends to the zero matrix. This happens if and only if $\rho(R) < 1$.

## 2. Proof of Necessity ($\rho(R) \ge 1 \implies$ Divergence)

Suppose $\rho(R) \ge 1$. Let $\lambda$ be an eigenvalue of $R$ such that $|\lambda| = \rho(R)$, and let $\mathbf{z} \neq \mathbf{0}$ be its corresponding eigenvector: $R\mathbf{z} = \lambda\mathbf{z}$.

If we choose the initial error $\mathbf{e}^{(0)}$ to be proportional to this eigenvector, $\mathbf{e}^{(0)} = \alpha\mathbf{z}$ ($\alpha \neq 0$), then the error at step $k$ is:

$$\mathbf{e}^{(k)} = R^k \mathbf{e}^{(0)} = R^k (\alpha\mathbf{z}) = \alpha (R^k \mathbf{z}) = \alpha \lambda^k \mathbf{z}$$

Taking the norm of the error vector:

$$\|\mathbf{e}^{(k)}\| = |\alpha| |\lambda|^k \|\mathbf{z}\| = |\alpha| (\rho(R))^k \|\mathbf{z}\|$$

Since $\rho(R) \ge 1$ and $\mathbf{z} \neq \mathbf{0}$, the norm $\|\mathbf{e}^{(k)}\|$ does not converge to zero. Therefore, convergence fails if $\rho(R) \ge 1$.

## 3. Proof of Sufficiency (Diagonalizable Case)

Assume $\rho(R) < 1$. If $R$ is **diagonalizable**, there exists a nonsingular matrix $S$ such that $R = S \Lambda S^{-1}$, where $\Lambda$ is the diagonal matrix containing the eigenvalues $\lambda_i$ of $R$.

The $k$-th power of $R$ is given by:

$$R^k = S \Lambda^k S^{-1}$$

Since $\rho(R) < 1$, we have $|\lambda_i| < 1$ for all eigenvalues. As $k \to \infty$, the power of each individual eigenvalue approaches zero, $\lim_{k\to\infty} \lambda_i^k = 0$. Consequently, the diagonal matrix $\Lambda^k$ converges to the zero matrix:

$$\lim_{k\to\infty} \Lambda^k = \mathbf{0}$$

Therefore, $R^k$ also converges to the zero matrix:

$$\lim_{k\to\infty} R^k = S \left( \lim_{k\to\infty} \Lambda^k \right) S^{-1} = S \mathbf{0} S^{-1} = \mathbf{0}$$

Since $\mathbf{e}^{(k)} = R^k \mathbf{e}^{(0)}$, the error converges to zero for any initial guess $\mathbf{e}^{(0)}$.

## 4. Proof of Sufficiency (General Case)

This proof segment, which demonstrates the convergence of splitting methods when the iteration matrix $G$ is not diagonalizable, relies on constructing a specialized matrix norm that is arbitrarily close to the spectral radius $\rho(G)$. This technique ensures that if $\rho(G) < 1$, the norm of $G$ is also strictly less than 1, guaranteeing the convergence of the sequence $G^k$.

### 1. Leveraging the Spectral Radius Property

A fundamental result in matrix analysis states that for any matrix $G \in \mathbb{C}^{n \times n}$ and any $\epsilon > 0$, there exists a subordinate matrix norm (or consistent norm) $\|\cdot\|_*$ such that $\|G\|_* < \rho(G) + \epsilon$.

Since we assume $\rho(G) < 1$, we can choose $\epsilon$ such that $\rho(G) + \epsilon < 1$. If we can construct a norm $\| \cdot \|_Z$ such that $\|G\|_Z < 1$, then because all induced matrix norms are **consistent** (satisfying $\|G^k\|_Z \le \|G\|_Z^k$), the powers of $G$ must converge to zero:

$$\lim_{k\to\infty} \|\mathbf{e}^{(k)}\|_Z = \lim_{k\to\infty} \|G^k \mathbf{e}^{(0)}\|_Z \le \lim_{k\to\infty} \|G\|_Z^k \|\mathbf{e}^{(0)}\|_Z = 0$$

### 2. Using the Schur Decomposition to Construct the Norm

To explicitly construct such a norm, we use the **Schur decomposition**. Every square matrix $G$ can be decomposed as:

$$G = Q T Q^H$$

where $Q$ is a **unitary matrix** ($Q^H Q = I$) and $T$ is an **upper triangular matrix**. The diagonal elements of $T$, denoted $t_{ii}$, are the eigenvalues $\lambda_i$ of $G$.

We define a diagonal scaling matrix $D$ with strictly positive entries:

$$D = \text{diag}(\delta^{-1}, \delta^{-2}, \dots, \delta^{-n})$$

where $\delta$ is a small positive parameter to be determined later.

We now analyze the similarity transformation of $T$ by $D$: $B = D T D^{-1}$.

The entries of $B$ are found to be:

$$[D T D^{-1}]_{ij} = t_{ij} \delta^{j-i}, \quad j \ge i$$

and 0 otherwise.

We evaluate the matrix $\infty$-norm of $B$, $\|B\|_{\infty}$, which is the maximum absolute row sum:

$$\| D T D^{-1} \|_{\infty} \le \max_i \{ |t_{ii}| + n \max_{j>i} \{ |t_{ij}| \delta^{j-i} \} \}$$

Since $\rho(G) < 1$, all diagonal entries $|t_{ii}| = |\lambda_i| \le \rho(G)$. By choosing $\delta$ sufficiently small, we can make the off-diagonal terms (where $j > i$ and $\delta^{j-i}$ is positive) arbitrarily small.

Specifically, since $\rho(G) < 1$, we can choose $\delta$ small enough to ensure that:

$$\| D T D^{-1} \|_{\infty} \le \rho(G) + \epsilon < 1$$

### 3. Defining the Norm and Concluding Convergence

Define the change of basis matrix $Z = D Q^H$. Since $D$ and $Q$ are nonsingular, $Z$ is nonsingular. We use $Z$ to define a vector norm:

$$\|\mathbf{x}\|_Z = \|Z\mathbf{x}\|_{\infty}$$

The corresponding induced operator norm for $G$ is given by $\|G\|_Z = \|Z G Z^{-1}\|_{\infty}$.

Substituting the Schur form of $G$ and the definition of $Z$ into the norm expression:

$$\| G \|_Z = \| (D Q^H) (Q T Q^H) (D Q^H)^{-1} \|_{\infty} = \| D T D^{-1} \|_{\infty}$$

We established in the previous step that, by choosing $\delta$ appropriately, $\| D T D^{-1} \|_{\infty} < 1$.

Therefore, we have constructed an induced norm $\| \cdot \|_Z$ such that:

$$\| G \|_Z < 1$$

Since $\| G^k \|_Z \le \| G \|_Z^k$, the error term $\mathbf{e}^{(k)} = G^k \mathbf{e}^{(0)}$ converges to zero as $k \to \infty$, thus proving that the splitting method converges for the general case if $\rho(G) < 1$.
```

## Conceptual Analogy

The iteration matrix $R$ acts like a filter applied repeatedly to the error $e^{(0)}$. For the iteration to converge, $R$ must shrink the error in all possible directions. The eigenvalues of $R$ represent the scaling factors applied in specific invariant directions. If the largest scaling factor (the spectral radius $\rho(R)$) is less than 1, every component of the error, regardless of how complex the geometry of the matrix $R$ is (diagonalizable or not, analogous to a perfectly spherical shrink or a distorted, complicated squeeze), is guaranteed to decay exponentially to zero over time.