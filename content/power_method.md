# The Power Method

Our first and most fundamental iterative algorithm is the **power method** (or *power iteration*). It is a simple but surprisingly effective method for finding the single, strictly dominant eigenvalue—the one with the largest magnitude—and its corresponding eigenvector.

While we will eventually move to more sophisticated methods, the power method provides the core intuition for many advanced algorithms, including the celebrated $QR$ iteration.

## Intuition: The Dominant Direction

Imagine the matrix $A$ as a linear transformation. When it acts on a vector, it stretches or shrinks it. The eigenvectors are the special directions that are *not* rotated, only scaled.

Now, consider what happens when we apply this transformation repeatedly to an arbitrary starting vector $\boldsymbol{q}$. We can think of $\boldsymbol{q}$ as a mixture of all the eigenvector components.

$$
\boldsymbol{q} = c_1 \boldsymbol{x}_1 + c_2 \boldsymbol{x}_2 + \cdots + c_n \boldsymbol{x}_n
$$

Each time we apply $A$, each component $c_i \boldsymbol{x}_i$ is scaled by its eigenvalue $\lambda_i$.

$$
A \boldsymbol{q} = c_1 \lambda_1 \boldsymbol{x}_1 + c_2 \lambda_2 \boldsymbol{x}_2 + \cdots + c_n \lambda_n \boldsymbol{x}_n
$$

If we apply $A$ $k$ times, we get:

$$
A^k \boldsymbol{q} = c_1 \lambda_1^k \boldsymbol{x}_1 + c_2 \lambda_2^k \boldsymbol{x}_2 + \cdots + c_n \lambda_n^k \boldsymbol{x}_n
$$

If one eigenvalue, say $\lambda_1$, is larger in magnitude than all others ($|\lambda_1| > |\lambda_2| \ge \dots$), its $k$-th power $\lambda_1^k$ will grow much faster than all other terms. After many iterations, the first term will dominate the sum, and all other components will become negligible in comparison.

$$
A^k \boldsymbol{q} \approx c_1 \lambda_1^k \boldsymbol{x}_1
$$

This means the resulting vector $A^k \boldsymbol{q}$ will be almost perfectly aligned with the dominant eigenvector $\boldsymbol{x}_1$.

## Derivation

Let's formalize this. We assume $A \in \mathbb{C}^{n \times n}$ is diagonalizable and has a single, strictly dominant eigenvalue $\lambda_1$:

$$
|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots \ge |\lambda_n|
$$

Let $\boldsymbol{x}_1, \dots, \boldsymbol{x}_n$ be the corresponding eigenvectors. We pick a starting vector $\boldsymbol{q}_0$ which we can express as a linear combination of these eigenvectors:

$$
\boldsymbol{q}_0 = \sum_{i=1}^n c_i \boldsymbol{x}_i
$$

We must assume that $c_1 \neq 0$. That is, our starting vector must have some component in the direction of the dominant eigenvector. In practice, a randomly chosen $\boldsymbol{q}_0$ makes this overwhelmingly likely.

Now, let's examine the vector $\boldsymbol{z}_k = A^k \boldsymbol{q}_0$:

$$
\boldsymbol{z}_k = A^k \left( \sum_{i=1}^n c_i \boldsymbol{x}_i \right) = \sum_{i=1}^n c_i A^k \boldsymbol{x}_i = \sum_{i=1}^n c_i \lambda_i^k \boldsymbol{x}_i
$$

We factor out the dominant term $\lambda_1^k$:

$$
\boldsymbol{z}_k = \lambda_1^k \left( c_1 \boldsymbol{x}_1 + \sum_{i=2}^n c_i \left(\frac{\lambda_i}{\lambda_1}\right)^k \boldsymbol{x}_i \right)
$$

Because of our dominance assumption, $|\lambda_i / \lambda_1| < 1$ for all $i \ge 2$. Therefore, as $k \to \infty$, the summation term decays exponentially:

$$
\lim_{k \to \infty} \left[ \sum_{i=2}^n c_i \left(\frac{\lambda_i}{\lambda_1}\right)^k \boldsymbol{x}_i \right] = \boldsymbol{0}
$$

This leaves us with:

$$
\boldsymbol{z}_k = A^k \boldsymbol{q}_0 \approx c_1 \lambda_1^k \boldsymbol{x}_1
$$

This shows that the sequence of vectors $\boldsymbol{z}_k$ converges to the subspace spanned by the dominant eigenvector, $\text{span}\{\boldsymbol{x}_1\}$.

## The Algorithm: Power Iteration

In practice, computing $A^k \boldsymbol{q}_0$ directly is numerically unstable. If $|\lambda_1| > 1$, the vector $\boldsymbol{z}_k$ will overflow, and if $|\lambda_1| < 1$, it will underflow to zero.

The solution is to **normalize** the vector at every single iteration. This keeps the vector's norm at $1$ and preserves only the directional information, which is what we care about. This leads to the **Power Iteration** algorithm.

**Algorithm: Power Iteration**

1.  Start with a random unit vector $\boldsymbol{q}$ (i.e., $\|\boldsymbol{q}\|_2 = 1$).
2.  Loop until convergence (e.g., until $\lambda$ stops changing):
      1. $\boldsymbol{z} = A \boldsymbol{q}$
      2. $\lambda = \boldsymbol{q}^H \boldsymbol{z}$ (Estimate eigenvalue using the **Rayleigh quotient**)
      3. $\boldsymbol{q} = \boldsymbol{z} / \|\boldsymbol{z}\|_2$ (Normalize for the next iteration)
3.  Return $\lambda$ (the dominant eigenvalue) and $\boldsymbol{q}$ (the dominant eigenvector).

## Convergence Analysis

* **Eigenvector:** The convergence of the vector $\boldsymbol{q}_k$ is slightly subtle. From our derivation, the normalized vector $\boldsymbol{q}_k = \boldsymbol{z}_k / \|\boldsymbol{z}_k\|$ behaves like:

    $$
    \boldsymbol{q}_k \approx \frac{c_1 \lambda_1^k \boldsymbol{x}_1}{\|c_1 \lambda_1^k \boldsymbol{x}_1\|} = \left( \frac{\lambda_1}{|\lambda_1|} \right)^k \left( \frac{c_1}{|c_1|} \right) \left( \frac{\boldsymbol{x}_1}{\|\boldsymbol{x}_1\|} \right)
    $$

    If $\lambda_1$ is real and positive, $\boldsymbol{q}_k$ converges to $\pm \boldsymbol{x}_1$. However, if $\lambda_1 = \rho e^{i\theta}$ is complex, the phase factor $(e^{i\theta})^k$ remains, and the vector $\boldsymbol{q}_k$ "spins" in the complex plane.

   * In all cases, the *subspace* $\text{span}\{\boldsymbol{q}_k\}$ converges to $\text{span}\{\boldsymbol{x}_1\}$.
   * The *rate of convergence* is linear and depends on the ratio $|\lambda_2 / \lambda_1|$. If this ratio is close to 1 (i.e., the dominant eigenvalue is not very dominant), convergence can be very slow.

* **Eigenvalue:** The sequence of estimates $\lambda_k = \boldsymbol{q}_{k}^H A \boldsymbol{q}_{k}$ converges to the dominant eigenvalue $\lambda_1$. We prove this result below.

````{prf:proof}
We prove that the sequence of eigenvalue estimates $\lambda_k = \boldsymbol{q}_k^H A \boldsymbol{q}_k$ (the Rayleigh quotient) converges to the dominant eigenvalue $\lambda_1$.

**Recall Eigenvector Convergence:**

From the analysis of the eigenvector, we know that as $k \to \infty$, the vector $\boldsymbol{q}_k$ becomes almost perfectly aligned with the dominant eigenvector $\boldsymbol{x}_1$.

Let's assume $\boldsymbol{x}_1$ is normalized such that $\|\boldsymbol{x}_1\|_2 = 1$. As shown in the derivation, the sequence $\boldsymbol{q}_k$ can be written as:

$$
\boldsymbol{q}_k = \phi_k \boldsymbol{x}_1 + \boldsymbol{\epsilon}_k
$$

where:

* $\phi_k = \left( \frac{\lambda_1}{|\lambda_1|} \right)^k \left( \frac{c_1}{|c_1|} \right)$ is a complex phase factor. Note that $|\phi_k| = 1$, which means $\phi_k^H \phi_k = 1$.
* $\boldsymbol{\epsilon}_k$ is the error vector, which contains all the non-dominant eigenvector components. From our derivation, we know that $\boldsymbol{\epsilon}_k \to \boldsymbol{0}$ as $k \to \infty$.

**Substitute into the Rayleigh Quotient:**

Now, let's examine the eigenvalue estimate $\lambda_k$ using this expression for $\boldsymbol{q}_k$:

$$
\lambda_k = \boldsymbol{q}_k^H A \boldsymbol{q}_k = (\phi_k \boldsymbol{x}_1 + \boldsymbol{\epsilon}_k)^H A (\phi_k \boldsymbol{x}_1 + \boldsymbol{\epsilon}_k)
$$

Since $A(\phi_k \boldsymbol{x}_1) = \phi_k (A \boldsymbol{x}_1) = \phi_k \lambda_1 \boldsymbol{x}_1$, we can write:

$$
\lambda_k = (\phi_k^H \boldsymbol{x}_1^H + \boldsymbol{\epsilon}_k^H) (\phi_k \lambda_1 \boldsymbol{x}_1 + A \boldsymbol{\epsilon}_k)
$$

**Expand the Expression:**

We multiply the terms:

$$
\lambda_k = \underbrace{(\phi_k^H \boldsymbol{x}_1^H)(\phi_k \lambda_1 \boldsymbol{x}_1)}_{\text{Term 1}} + \underbrace{(\phi_k^H \boldsymbol{x}_1^H)(A \boldsymbol{\epsilon}_k)}_{\text{Term 2}} + \underbrace{(\boldsymbol{\epsilon}_k^H)(\phi_k \lambda_1 \boldsymbol{x}_1)}_{\text{Term 3}} + \underbrace{(\boldsymbol{\epsilon}_k^H A \boldsymbol{\epsilon}_k)}_{\text{Term 4}}
$$

**Analyze Each Term as $k \to \infty$:**

* **Term 1:** We can regroup the scalars:

    $$
    (\phi_k^H \phi_k) \lambda_1 (\boldsymbol{x}_1^H \boldsymbol{x}_1)
    $$

    Since $\phi_k^H \phi_k = 1$ (as $|\phi_k|=1$) and $\boldsymbol{x}_1^H \boldsymbol{x}_1 = 1$ (as $\boldsymbol{x}_1$ is normalized), this term simplifies to:

    $$
    \text{Term 1} = (1) \cdot \lambda_1 \cdot (1) = \lambda_1
    $$

* **Term 2:** As $k \to \infty$, we know $\boldsymbol{\epsilon}_k \to \boldsymbol{0}$. The entire term, which is a scalar inner product, must therefore go to 0.

* **Term 3:** Similarly, as $k \to \infty$, this inner product also goes to 0.

* **Term 4:** This is a quadratic form $\boldsymbol{\epsilon}_k^H A \boldsymbol{\epsilon}_k$. This term goes to 0 even faster.

By taking the limit of the expanded expression, we get:

$$
\begin{aligned}
\lim_{k\to\infty}\lambda_k
&= \lim_{k\to\infty}\big(\text{Term 1} + \text{Term 2} + \text{Term 3} + \text{Term 4}\big) \\
&= \lambda_1 + 0 + 0 + 0 \\
&= \lambda_1
\end{aligned}
$$

Therefore,

$$
\lim_{k \to \infty} \lambda_k = \lambda_1
$$

This proves that the sequence of Rayleigh quotients converges to the dominant eigenvalue $\lambda_1$.
````

## Limitations

The power method is fundamental, but it has significant limitations:

1.  It only finds the **single largest eigenvalue**.
2.  It fails if there is no *single* dominant eigenvalue (e.g., a complex conjugate pair $|\lambda_1| = |\lambda_2|$ or real eigenvalues $\lambda_1 = -\lambda_2$).
3.  Convergence can be unacceptably slow if $|\lambda_1| \approx |\lambda_2|$.

Despite these limitations, this core idea of "repeated application and normalization" is the engine inside the much more robust algorithms we will study next.