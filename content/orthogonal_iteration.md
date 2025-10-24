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

It is important to note that $Q$ in the Schur decomposition is not unique. If $A = Q T Q^H$, we can define a new $Q' = Q D$, where $D$ is any diagonal matrix with entries $e^{i \theta_j}$ ($\theta_j \in \mathbb R$, $|e^{i \theta_j}|=1$). $Q'$ is still unitary, and 

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

So, $\text{span}\{\boldsymbol{q}_{k,1}\} \to \text{span}\{\boldsymbol{v}_1\}$ as $k \to \infty$. The base case $P(1)$ holds.

### Inductive Step: Assume $P(i)$ holds, prove $P(i+1)$

**Assumption $P(i)$:** We assume that $\mathcal{S}_{k,i} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,i}\}$ converges to $\mathcal{V}_i = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_i\}$.

**Goal:** We must show that $P(i+1)$ holds, i.e., $\mathcal{S}_{k,i+1} = \text{span}\{\boldsymbol{q}_{k,1}, \dots, \boldsymbol{q}_{k,i+1}\}$ converges to $\mathcal{V}_{i+1} = \text{span}\{\boldsymbol{v}_1, \dots, \boldsymbol{v}_{i+1}\}$.

Let's analyze the $(i+1)$-th column of the iteration, $\boldsymbol{q}_{k+1, i+1}$. The $QR$ factorization $Q_{k+1} R_{k+1} = Z$ is equivalent to the Gram-Schmidt process. The vector $\boldsymbol{q}_{k+1, i+1}$ is computed by taking $\boldsymbol{z}_{i+1} = A \boldsymbol{q}_{k, i+1}$ and orthogonalizing it against the *preceding* vectors $\boldsymbol{q}_{k+1, 1}, \dots, \boldsymbol{q}_{k+1, i}$.

$$
\begin{aligned}
\boldsymbol{w} &= A \boldsymbol{q}_{k, i+1} - \sum_{j=1}^{i} (\boldsymbol{q}_{k+1, j}^H (A \boldsymbol{q}_{k, i+1})) \boldsymbol{q}_{k+1, j} \\
\boldsymbol{q}_{k+1, i+1} &= \frac{\boldsymbol{w}}{\|\boldsymbol{w}\|_2}
\end{aligned}
$$

The vector $\boldsymbol{w}$ is the component of $A \boldsymbol{q}_{k, i+1}$ that is orthogonal to $\mathcal{S}_{k+1, i}$.

1. Let $P^{(i)}$ be the true projection onto the orthogonal complement of $\mathcal{V}_i$, i.e., $P^{(i)} = I - \sum_{j=1}^i \boldsymbol{v}_j \boldsymbol{v}_j^H$. Let $P_k^{(i)}$ be the projection onto the orthogonal complement of $\mathcal{S}_{k,i}$.

2. From our inductive assumption $P(i)$, we know $\mathcal{S}_{k,i} \to \mathcal{V}_i$. This means the projection $P_k^{(i)} \to P^{(i)}$ as $k \to \infty$.

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


## Second Proof of Convergence


We provide another convergence proof that is not based on the use of orthogonal projections and the idea of deflation. This proof is more precise but also significantly more technical. It will give us a precise convergence rate for the method depending on the ratio $\lambda_{p+1} / \lambda_p$.

```{admonition} Summary of key ideas and results
:class: tip

- **Coordinates:** Work in a Schur basis $A=UTU^\ast$ and express the current iterate $U^\ast Q_k$ (where $Q_k$ is assumed to have $p$ columns) as a graph 

$$\mathrm{span}\!\begin{bmatrix} C_k \\ S_k \end{bmatrix} = \mathrm{span}\!\begin{bmatrix}I\\E_k\end{bmatrix}$$ 

with $E_k=S_k C_k^{-1}$.
- **One-step update:** The graph variable obeys

$$
E_{k+1}=T_{22}\,E_k\,(T_{11}+T_{12}E_k)^{-1}.
$$

- **Local convergence:** If the target and unwanted spectra are separated (formally, $\mathrm{sep}(T_{11},T_{22})>0$) and the start has nonzero overlap with the target subspace ($C_0$ invertible), then $E_k\to 0$ **linearly**; asymptotically the factor is about $\rho(T_{22})\,\rho(T_{11}^{-1})$.
- **Global rate with a modulus gap:** If $|\lambda_p|>|\lambda_{p+1}|$, then the error decays roughly like

$$
\|E_k\|\;\lesssim\; C\left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k,
$$

up to constants.
- **What converges:** The subspace $\mathrm{span}(Q_k)$ converges to the Schur subspace $\mathrm{span}(U_1)$ and the Ritz values of $B_k:=Q_k^\ast A Q_k$ converge to the eigenvalues in $T_{11}$.
- **How to measure error:** $\|E_k\|=\|\tan\Theta_k\|$, where $\Theta_k$ are the principal angles (see below for a definition of principal angles) between $\mathrm{span}(Q_k)$ and $\mathrm{span}(U_1)$; hence those angles go to $0$ at the same rate.
```

### Setup (ordered Schur form and the iteration)

Let $A\in\mathbb{C}^{n\times n}$. Orthogonal iteration with block size $p$ (i.e., $Q_k$ has $p$ columns) is

$$
Y_{k+1}=A\,Q_k,\qquad Q_{k+1}=\operatorname{orth}(Y_{k+1}),\quad Q_0^\ast Q_0=I_p.
$$

The notation $\operatorname{orth}(\cdot)$ means to compute an orthonormal basis for the range (e.g., via $QR$).

We will work with the (unitary) Schur form. Choose a unitary $U$ so that

$$
A=U\begin{bmatrix}T_{11}&T_{12}\\0&T_{22}\end{bmatrix}U^\ast=:U T U^\ast,
$$

where the $p\times p$ block $T_{11}$ holds the eigenvalues we want (e.g., the $p$ largest in modulus, or any isolated cluster after a reordering), and $T_{22}$ holds the rest. The convergence we wish to analyze is toward the Schur subspace $\mathcal{U}_1=\operatorname{span}(U_1)$ (columns of $U$ corresponding to $T_{11}$).

Write the iterate in Schur coordinates as $Z_k:=U^\ast Q_k$ and block it

$$
Z_k=\begin{bmatrix}C_k\\ S_k\end{bmatrix},\qquad Z_k^\ast Z_k=I_p.
$$

Subspace convergence can be studied through $Z_k$.

### The graph map and the core recurrence

As long as $C_k$ is nonsingular (true once the iterate has nonzero overlap with $\mathcal{U}_1$), define the graph variable $E_k:=S_k C_k^{-1}$ so that 

$$
\operatorname{span}(Z_k)=\operatorname{span}\!\begin{bmatrix}I\\ E_k\end{bmatrix}.
$$

The name graph variable comes from functional analysis and corresponds to viewing the subspace as the graph of a linear operator from the "top" to the "bottom" space.

One step of the (unnormalized) iteration in Schur coordinates is

$$
U^\ast Y_{k+1}=T Z_k=\begin{bmatrix}T_{11}C_k+T_{12}S_k\\ T_{22}S_k\end{bmatrix}.
$$

Since column spaces are what matter, the next subspace can be represented as a graph too:

$$
E_{k+1} = T_{22} S_k \bigl(T_{11} C_k + T_{12} S_k\bigr)^{-1} =
T_{22}\,E_k\,\bigl(T_{11}+T_{12}E_k\bigr)^{-1}
$$

using $S_k=E_k C_k$.

This “graph transform” is the standard invariant-subspace map for block-upper-triangular $T$; $E=0$ is its fixed point (the exact invariant subspace $\mathcal{U}_1$). Linearizing at $E=0$ gives

$$
E_{k+1}=T_{22}\,E_k\,T_{11}^{-1}+ \mathcal{O}(\|E_k\|^2).
$$

Hence, locally the contraction factor is governed by $T_{22}\,(\cdot)\,T_{11}^{-1}$.

### Local linear convergence

From the graph recurrence and a Neumann-series bound,

$$
\|E_{k+1}\|\;\le\;\|T_{22}\|\,\|T_{11}^{-1}\|\;\frac{\|E_k\|}{\,1-\|T_{11}^{-1}\|\,\|T_{12}\|\,\|E_k\|\,}.
$$

Therefore, once $\|E_k\|$ is small enough (in particular, $\|T_{11}^{-1}\|\,\|T_{12}\|\,\|E_k\|<1$), you get linear contraction with asymptotic factor $\|T_{22}\|\,\|T_{11}^{-1}\|$. This is the standard local behavior around $E=0$ in the non-Hermitian case.

### Global convergence and rates (dominant-modulus case)

Assume the eigenvalues are ordered by modulus

$$
|\lambda_1| > \cdots > |\lambda_p| \;>\; |\lambda_{p+1}| > \cdots > |\lambda_n|,
$$

and the initial subspace has a nonzero component in $\mathcal{U}_1$. Then the orthogonal-iteration subspaces $\operatorname{span}(Q_k)$ converge to $\mathcal{U}_1$ (the Schur subspace for $\lambda_1,\ldots,\lambda_p$). Moreover, the principal-angle (or projector) errors decay essentially like powers of the modulus gap; asymptotically, each column converges at a rate $|\lambda_{i+1}/\lambda_i|$ (for $i=1,\dots,p$):

$$
\sin\angle\bigl(\mathcal{U}_1,\operatorname{span}(Q_k)\bigr)\;\lesssim\; \kappa \, \Bigl(\tfrac{|\lambda_{p+1}|}{|\lambda_p|}\Bigr)^{k}.
$$

Here $\kappa$ reflects conditioning of the eigenbasis (departure from normality). This is the rate one gets from the linearized map $E_{k+1}\approx T_{22}E_k T_{11}^{-1}$.

```{prf:theorem} Informal theorem
:label: thm:orth-it-nh

If $\Lambda(T_{11})$ and $\Lambda(T_{22})$ are disjoint and $|\lambda_p|>|\lambda_{p+1}|$, then for generic $Q_0$ the iterates $\operatorname{span}(Q_k)$ converge to $\mathcal{U}_1$. Asymptotically, principal angles decay at a linear rate determined by $|\lambda_{p+1}/\lambda_p|$, and locally the graph error satisfies $\|E_{k+1}\|\le \|T_{22}\|\,\|T_{11}^{-1}\|\,\|E_k\|+o(\|E_k\|)$.
```

```{prf:theorem} Convergence to a Schur invariant subspace
:label: thm:orth-it-nh-rigorous

(a) (Local convergence under spectral separation)  
Assume the **separation**

$$
\mathrm{sep}(T_{11},T_{22}) \;:=\; \min_{X\neq 0} \frac{\|T_{11}X - X T_{22}\|}{\|X\|} \;>\; 0.
$$

Then there exists $\varepsilon>0$ such that if $\|E_0\|<\varepsilon$ the orthogonal iteration is well-defined for all $k$ (meaning $C_k$ is invertible so $E_k=S_k C_k^{-1}$ exists, and $T_{11}+T_{12}E_k$ is invertible so the update for $E_{k+1}$ is valid) and the graph iterates satisfy the Riccati map

$$
E_{k+1} = T_{22}\,E_k\,(T_{11}+T_{12}E_k)^{-1}.
$$

Moreover, there are constants $c$, $\delta>0$ so that for all $\|E_k\|<\delta$,

$$
\|E_{k+1}\|\;\le\;\|T_{22}\|\,\|T_{11}^{-1}\|\,\|E_k\| \;+\; c\,\|E_k\|^2,
$$

hence $E_k\to 0$ **linearly**. The asymptotic rate obeys

$$
\limsup_{k\to\infty}\frac{\|E_{k+1}\|}{\|E_k\|} \;\le\;
\rho(T_{22})\,\rho(T_{11}^{-1}),
$$

where $\rho(\cdot)$ is the spectral radius. Writing $\Theta_k$ for the principal-angle matrix between $\mathrm{span}(Q_k)$ and $\mathrm{span}(U_1)$, we have $\|\tan\Theta_k\| = \|E_k\|$ and $\|\sin\Theta_k\|\le \|E_k\|$, so the subspaces converge.

(b) (Global rate under a modulus gap)  
Assume the eigenvalues are ordered by modulus with a gap

$$
|\lambda_1| > \cdots > |\lambda_p| \;>\; |\lambda_{p+1}| > \cdots >|\lambda_n|.
$$

If $C_0$ is nonsingular, then the subspaces $\mathrm{span}(Q_k)$ converge to the Schur subspace $\mathrm{span}(U_1)$. There exist constants $c>0$ such that, for all sufficiently large $k$,

$$
\|E_k\| \;\le\; c\,\left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k.
$$

Consequently, the largest principal angle satisfies

$$
\sin\theta_{\max}(\mathrm{span}(Q_k),\mathrm{span}(U_1))
\;\le\; c\,\left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k.
$$

Finally, the Ritz matrix $B_k:=Q_k^\ast A Q_k$ has eigenvalues that converge to the eigenvalues of $T_{11}$.
```

### Practical notes

- **Non-normality can slow convergence.** Poorly conditioned eigenvectors (large departure from normality) can introduce substantial polynomial factors before the geometric rate dominates. It is one reason subspace iteration targets Schur vectors and/or uses accelerations.
- **Acceleration:** shift-and-invert $(A-\sigma I)^{-1}$, polynomial filters (Chebyshev), and locking/deflation are standard to isolate interior clusters and improve gaps; these plug directly into the subspace iteration framework.
- **Separation and Sylvester equations:** the spectral separation $\mathrm{sep}(T_{11},T_{22})>0$ ensures a well-conditioned invariant-subspace problem and appears in error/conditioning bounds (e.g., Schur–Parlett/sign function analyses).


### Definition of Principal Angles

We define the concept of principal angles between two subspaces, which were used in the convergence analysis above.

Let $\mathcal U$ and $\mathcal V$ be subspaces of $\mathbb{C}^n$ (or $\mathbb{R}^n$) with dimensions $p$ and $q$, where $p \le q$.

There exist orthonormal bases

$$
\{u_1, \dots, u_p\} \text{ for } \mathcal U,
\qquad
\{v_1, \dots, v_q\} \text{ for } \mathcal V,
$$

and real numbers

$$
0 \le \theta_1 \le \theta_2 \le \cdots \le \theta_p \le \tfrac{\pi}{2}
$$

such that

$$
\langle u_i, v_j \rangle =
\begin{cases}
\cos(\theta_i), & \text{if } i = j, \\
0, & \text{if } i \ne j.
\end{cases}
$$

The numbers $\theta_1, \dots, \theta_p$ are called the **principal angles** (or **canonical angles**) between $\mathcal U$ and $\mathcal V$.

Alternatively, if $U \in \mathbb{C}^{n \times p}$ and $V \in \mathbb{C}^{n \times q}$ have orthonormal columns spanning $\mathcal U$ and $\mathcal V$, then the cosines of the principal angles are the singular values of $U^* V$:

$$
\cos(\theta_i) = \sigma_i(U^*V), \quad i = 1, \dots, p.
$$

```{prf:lemma} Relationship with projector norms
Let $P_{\mathcal U}$ and $P_{\mathcal V}$ denote the orthogonal projectors onto $\mathcal U$ and $\mathcal V$, respectively, and let $\Theta = \operatorname{diag}(\theta_1, \dots, \theta_p)$ be the diagonal matrix of principal angles.

Then:

$$
\|\sin \Theta\| = \| P_{\mathcal U^\perp} P_{\mathcal V} \|
= \| P_{\mathcal U} - P_{\mathcal V} \|,
$$

and

$$
\|\tan \Theta\| = \| P_{\mathcal U^\perp} P_{\mathcal V} (P_{\mathcal U} P_{\mathcal V})^{-1} \|
$$

whenever the intersection $\mathcal U \cap \mathcal V$ is trivial.

In particular, $\sin \Theta$ and $\tan \Theta$ provide quantitative measures of the distance or “gap” between the two subspaces.
```

**Remarks**

* $\theta_1=0$ if and only if $\mathcal U \cap \mathcal V \neq \{0\}$.
* The largest principal angle $\theta_p$ provides a measure of the worst‐alignment of $\mathcal U$ with $\mathcal V$.
* In many numerical analyses (e.g., subspace iteration, perturbation bounds), one works with $\sin\theta_i$ or $\tan\theta_i$ rather than $\theta_i$ itself.  ￼
* When $p=q=1$ (lines), this notion reduces to the usual angle between two lines.

