# Convergence of the Arnoldi Process

The [previous section](arnoldi_iteration.md) constructed the orthonormal basis $Q_k$ and its projected matrix $H_k$. To understand the resulting eigenvalue approximations, we need to connect the directions in the subspace to the Ritz pairs extracted from it.

Polynomials provide that connection. Every vector in a Krylov subspace has the form $p(A)q_1$. We will show how the Ritz values determine such polynomials, how these polynomials can isolate a wanted eigenvector, and how the Arnoldi residual measures the accuracy of the resulting pair.

We work in exact arithmetic with $A\in\mathbb C^{n\times n}$ and $\|q_1\|_2=1$. For real matrices and vectors, conjugate transpose reduces to transpose.

Recall the Arnoldi relation

$$
AQ_k=Q_kH_k+f_ke_k^T,
\qquad Q_k^Hf_k=0.
$$

Here $e_k$ selects the last column, and $f_k$ is the part of $Aq_k$ perpendicular to the current subspace:

$$
f_k=Aq_k-Q_kH_ke_k.
$$

If this remainder is nonzero, $f_k=h_{k+1,k}q_{k+1}$, where $h_{k+1,k}=\|f_k\|_2$ and $q_{k+1}$ is a unit vector. If $f_k=0$, the subspace is invariant.

(sec:arnoldi-ritz-residuals)=
## Extracting Ritz Pairs and Checking Their Residuals

Apply the [Rayleigh–Ritz construction](krylov_subspaces.md) to the matrix produced by Arnoldi:

$$
H_ky=\theta y,\qquad \|y\|_2=1.
$$

The associated Ritz vector is $u=Q_ky$, with $\|u\|_2=1$. A single eigenpair of $H_k$ suffices; $H_k$ need not be diagonalizable.

Define the residual $r=Au-\theta u$ to measure how closely this pair satisfies the original eigenvector equation. Multiplying the Arnoldi relation by $y$ gives

$$
\begin{aligned}
Au
&=AQ_ky\\
&=Q_kH_ky+f_k(e_k^Ty)\\
&=\theta Q_ky+f_k(e_k^Ty)\\
&=\theta u+f_k(e_k^Ty).
\end{aligned}
$$

Therefore,

$$
\boxed{
\begin{aligned}
r&=f_k(e_k^Ty),\\
\|r\|_2&=\|f_k\|_2\,|e_k^Ty|
=|h_{k+1,k}|\,|(y)_k|.
\end{aligned}
}
$$

Here $e_k^Ty=(y)_k$ is the last component of $y$, and $\|f_k\|_2=h_{k+1,k}$. This check requires only the remainder norm and one component of the small eigenvector. We can assess each Ritz pair without forming $u$ or performing another multiplication by $A$.

In Python, all the residual estimates can be computed together:

```python
# Use arnoldi_iteration from the preceding section.
Q, H, f = arnoldi_iteration(A, q1, k)
theta, Y = np.linalg.eig(H)

# Normalize the small eigenvectors before comparing residuals.
Y /= np.linalg.norm(Y, axis=0)
residual_estimates = np.linalg.norm(f) * np.abs(Y[-1, :])

# Form a particular Ritz vector only when it is needed.
j = 0
u = Q @ Y[:, j]
actual_residual = np.linalg.norm(A @ u - theta[j] * u)
```

The estimates use the exact-arithmetic identity and can differ from explicitly computed residuals because of rounding. Near the requested accuracy, an explicit residual check is useful. Residuals should be judged relative to the scale of the problem.

This gives us a practical accuracy check. To explain why increasing $k$ can improve the approximations, we next examine the polynomial combinations available in the Krylov subspace. We return to the eigenvalue error guaranteed by a small residual after that analysis.

## Polynomials Describe the Available Approximations

The basis result proved on the preceding page, {prf:ref}`thm:arnoldi_krylov_basis`, lets us analyze Arnoldi through the generators of its Krylov subspace. A linear combination of those vectors can be written as

$$
\begin{aligned}
v&=a_0q_1+a_1Aq_1+\cdots+a_{k-1}A^{k-1}q_1\\
 &=p(A)q_1,
\qquad p(z)=a_0+a_1z+\cdots+a_{k-1}z^{k-1}.
\end{aligned}
$$

Consequently,

$$
\boxed{
\mathcal K_k(A,q_1)
=\{p(A)q_1:\deg p\le k-1\}.
}
$$

The power method uses one particular polynomial, $p(z)=z^{k-1}$. Arnoldi makes the whole space of these polynomial combinations available to the projected eigenvalue problem.

### The Cayley–Hamilton Theorem

We want to relate polynomials of $A$ to the eigenvalues computed from $H_k$. Cayley–Hamilton gives the starting point: a matrix's characteristic polynomial annihilates that matrix. We will apply this fact to $H_k$ to derive Arnoldi's polynomial optimality property.

:::{prf:definition} Characteristic Polynomial
:label: def:characteristic_polynomial

The characteristic polynomial of $A$ is the monic polynomial

$$
p_A(z)=\det(zI-A).
$$

Here *monic* means that its leading coefficient is one.
:::

:::{prf:theorem} Cayley–Hamilton Theorem
:label: thm:cayley_hamilton

Every square matrix satisfies

$$
p_A(A)=0.
$$
:::

:::{prf:proof}
Use a complex Schur decomposition $A=QTQ^H$, where $T$ is upper triangular with diagonal $\lambda_1,\ldots,\lambda_n$. Then

$$
p_A(z)=\prod_{j=1}^n(z-\lambda_j).
$$

Let $E_j=\operatorname{span}\{e_1,\ldots,e_j\}$, with $E_0=\{0\}$. Upper triangularity gives

$$
(T-\lambda_jI)E_j\subseteq E_{j-1}.
$$

Indeed, multiplication cannot create a coordinate beyond $j$, and the $j$th coordinate becomes zero. Applying the factors from right to left therefore sends $E_n$ successively into $E_{n-1},\ldots,E_0$:

$$
p_A(T)=(T-\lambda_1I)\cdots(T-\lambda_nI)=0.
$$

Hence $p_A(A)=Qp_A(T)Q^H=0$. This proof includes matrices that are not diagonalizable.
:::

Cayley–Hamilton guarantees an annihilating polynomial of degree $n$. A lower-degree polynomial may already make $p(A)=0$, which explains why Arnoldi can reach an invariant subspace before exploring all $n$ directions.

:::{prf:definition} Minimal Polynomial
:label: def:minimal_polynomial

The minimal polynomial $m_A$ is the monic polynomial of smallest degree satisfying $m_A(A)=0$.
:::

:::{prf:proposition} Minimal Polynomial of a Diagonalizable Matrix
:label: prop:minimal_polynomial

If $A$ is diagonalizable and its distinct eigenvalues are $\nu_1,\ldots,\nu_d$, then

$$
m_A(z)=\prod_{i=1}^d(z-\nu_i).
$$
:::

:::{prf:proof}
For $A=X\Lambda X^{-1}$, the identity $p(A)=Xp(\Lambda)X^{-1}$ shows that $p(A)=0$ exactly when $p$ vanishes at every distinct eigenvalue.
:::

If $m_A$ has degree $d$, the identity $m_A(A)q_1=0$ expresses $A^dq_1$ using earlier powers. Arnoldi must therefore reach an invariant subspace in at most $d$ steps in exact arithmetic. It may stop earlier if a lower-degree polynomial annihilates the particular starting vector.

## Arnoldi's Polynomial Optimality Property

Exact termination corresponds to a polynomial that annihilates the starting vector. Before then, we can seek the smallest possible $\|p(A)q_1\|_2$ among **monic degree-$k$ polynomials**. The characteristic polynomial of $H_k$ attains this minimum. Its roots are the Ritz values, so the result also gives a way to understand their locations.

:::{prf:theorem} Arnoldi Optimality
:label: thm:arnoldi_optimality

Suppose Arnoldi has constructed $Q_k$ and $H_k$, with no breakdown before step $k$. Let

$$
p_k(z)=\det(zI-H_k).
$$

Then $p_k$ minimizes $\|p(A)q_1\|_2$ among all monic polynomials of degree $k$.
:::

:::{prf:proof}
Use the Arnoldi relation recalled above, where $f_k$ is the remainder outside the current subspace. Multiplication by an upper Hessenberg matrix can extend the nonzero entries of a vector downward by at most one position. Starting from $e_1$, the vector $H_k^me_1$ therefore has zeros beyond position $m+1$. In particular,

$$
e_k^TH_k^me_1=0,\qquad 0\le m\le k-2.
$$

Starting from $q_1=Q_ke_1$, suppose $A^mq_1=Q_kH_k^me_1$. For these values of $m$, the next multiplication gives

$$
\begin{aligned}
A^{m+1}q_1
&=(AQ_k)H_k^me_1\\
&=Q_kH_k^{m+1}e_1+f_k\underbrace{(e_k^TH_k^me_1)}_{0}\\
&=Q_kH_k^{m+1}e_1.
\end{aligned}
$$

This proves by induction that

$$
A^mq_1=Q_kH_k^me_1,\qquad 0\le m<k.
$$

At $m=k-1$, the last coordinate of $H_k^{k-1}e_1$ can be nonzero. The same multiplication now retains the remainder term:

$$
A^kq_1=Q_kH_k^ke_1+\gamma_kf_k,
\qquad
\gamma_k=e_k^TH_k^{k-1}e_1.
$$

Since $p_k$ is monic and $p_k(H_k)=0$ by Cayley–Hamilton,

$$
p_k(A)q_1
=Q_kp_k(H_k)e_1+\gamma_kf_k
=\gamma_kf_k\perp\mathcal K_k.
$$

For any other monic degree-$k$ polynomial $p$, the difference $s=p-p_k$ has degree at most $k-1$, so $s(A)q_1\in\mathcal K_k$. Orthogonality now gives

$$
\begin{aligned}
\|p(A)q_1\|_2^2
&=\|p_k(A)q_1+s(A)q_1\|_2^2\\
&=\|p_k(A)q_1\|_2^2+\|s(A)q_1\|_2^2\\
&\ge\|p_k(A)q_1\|_2^2.
\end{aligned}
$$
:::

The theorem minimizes a vector norm. We next translate that result into information about the polynomial's roots, then remove one root factor to recover an individual Ritz-vector direction.

## Connection to Eigenvalue Approximation

### Why a Small Polynomial Norm Matters

Write the roots of $p_k$ as $\theta_1,\ldots,\theta_k$, counting multiplicity. These are the Ritz values, so

$$
p_k(z)=\prod_{\ell=1}^k(z-\theta_\ell).
$$

Suppose $A=X\Lambda X^{-1}$, with unit-norm eigenvectors as columns of $X$, and $q_1=Xc$. Then

$$
X^{-1}p_k(A)q_1
=\operatorname{diag}\bigl(p_k(\lambda_1),\ldots,p_k(\lambda_n)\bigr)c.
$$

Consequently, for each nonzero starting coefficient $c_i$,

$$
|p_k(\lambda_i)|
\le
\frac{\|X^{-1}\|_2}{|c_i|}\,\|p_k(A)q_1\|_2.
$$

A small polynomial norm therefore makes $p_k(\lambda_i)$ small when this eigenvector component is well represented and the eigenvector coordinates are not too sensitive. Since $p_k(\lambda_i)$ is a product of distances to Ritz values, at least one of those distances must be small. More precisely,

$$
\boxed{
\min_\ell|\lambda_i-\theta_\ell|
\le
\left(
\frac{\|X^{-1}\|_2}{|c_i|}
\|p_k(A)q_1\|_2
\right)^{1/k}.
}
$$

Indeed, if every distance is at least $\delta$, their product is at least $\delta^k$. This explains the relevance of polynomial optimality to eigenvalue approximation. The $k$th root can make this a weak bound, however, and a tiny $c_i$ or sensitive eigenvector basis weakens it further. It also does not say that every Ritz value is accurate. The individual residual checks remain the practical test.

### From Roots to Eigenvector Filters

The polynomial $p_k$ describes all $k$ Ritz values together. Removing one root factor gives the corresponding Ritz-vector direction and an exact expression for its residual.

:::{prf:theorem} Polynomial Representation and Residual of a Ritz Vector
:label: thm:arnoldi-polynomial-ritz-residual

Suppose Arnoldi is performed in exact arithmetic with $\|q_1\|_2=1$ and no breakdown before step $k$. Let $H_k$ be its projected Hessenberg matrix and $p_k(z)=\det(zI-H_k)$. Write $\theta_1,\ldots,\theta_k$ for the eigenvalues of $H_k$, counting multiplicity. For any $1\le\ell\le k$, define

$$
s_\ell(z)=\frac{p_k(z)}{z-\theta_\ell}
=\prod_{j\ne\ell}(z-\theta_j).
$$

Then $v_\ell=s_\ell(A)q_1$ is nonzero, and its unit normalization

$$
u_\ell=\frac{v_\ell}{\|v_\ell\|_2}
$$

is a Ritz vector associated with $\theta_\ell$. Its residual satisfies

$$
\boxed{
\|Au_\ell-\theta_\ell u_\ell\|_2
=
\frac{\|p_k(A)q_1\|_2}{\|s_\ell(A)q_1\|_2}.
}
$$
:::

:::{prf:proof}
The polynomial $s_\ell$ has degree $k-1$, so $v_\ell=s_\ell(A)q_1$ lies in $\mathcal K_k$. It is nonzero: otherwise a monic degree-$(k-1)$ polynomial would make the first $k$ Krylov generators linearly dependent, contradicting the absence of earlier breakdown. The factorization gives

$$
(A-\theta_\ell I)v_\ell=p_k(A)q_1.
$$

The right side is perpendicular to $\mathcal K_k$ by the optimality proof. Since $v_\ell$ lies in that subspace, projecting gives

$$
(H_k-\theta_\ell I)Q_k^Hv_\ell
=Q_k^Hp_k(A)q_1=0.
$$

Thus $v_\ell$ is a Ritz-vector direction for $\theta_\ell$. Dividing $(A-\theta_\ell I)v_\ell=p_k(A)q_1$ by $\|v_\ell\|_2$ and taking norms proves the residual formula.
:::

Thus $\|p_k(A)q_1\|_2$ must be small relative to $\|s_\ell(A)q_1\|_2$, the norm of the vector before normalization, for this Ritz pair to have a small residual.

To compare this direction with a wanted eigenvector $x_1$ of eigenvalue $\lambda_1$, assume $s_\ell(\lambda_1)\ne0$ and rescale the same polynomial as

$$
p(z)=\frac{s_\ell(z)}{s_\ell(\lambda_1)}.
$$

This normalization preserves the wanted eigenvector component because $p(\lambda_1)=1$. The other Ritz values remain zeros of $p$, so they can suppress unwanted components. Rescaling changes only the length and phase of $p(A)q_1$, not its direction. The next bound quantifies when this direction is close to $x_1$; it does not assume in advance that the Ritz pair is accurate.

These polynomials are tools for analysis. The algorithm computes the Ritz pairs from $H_k$ without knowing $\lambda_1$ or explicitly forming $p_k$ and $s_\ell$.

### How a Polynomial Isolates an Eigenvector

Keep the diagonalizable matrix $A=X\Lambda X^{-1}$ and its unit-norm eigenvectors. To analyze the normalized filter, expand the starting vector as

$$
q_1=\sum_{i=1}^n c_ix_i.
$$

Applying the filter gives

$$
p(A)q_1=\sum_{i=1}^n c_i p(\lambda_i)x_i.
$$

Since $p(\lambda_1)=1$, the wanted component $c_1x_1$ is unchanged. The other components are multiplied by $p(\lambda_i)$, so small values at the unwanted eigenvalues make $p(A)q_1/c_1$ close to $x_1$, provided $c_1\ne0$. The following bound controls this particular vector as well as the best approximation in the subspace. It applies to the Ritz-root filter and to any other polynomial with the same degree and normalization.

:::{prf:theorem} Polynomial Bound for Eigenvector Approximation
:label: thm:arnoldi-polynomial-approximation

Suppose $n\ge2$, $A=X\Lambda X^{-1}$, the columns of $X$ have unit norm, and $q_1=Xc$ with $c_1\ne0$. Let $X_2=[x_2,\ldots,x_n]$ and $c_{\mathrm{rest}}=(c_2,\ldots,c_n)^T$ denote the unwanted eigenvectors and their coefficients, respectively.

For any polynomial $p$ of degree at most $k-1$ with $p(\lambda_1)=1$, set $v=p(A)q_1/c_1$. Then

$$
\boxed{
\begin{aligned}
\|(I-Q_kQ_k^H)x_1\|_2
&\le\|x_1-v\|_2,\\
\|x_1-v\|_2
&\le
\frac{\|X_2\|_2\|c_{\mathrm{rest}}\|_2}{|c_1|}
\max_{i>1}|p(\lambda_i)|.
\end{aligned}
}
$$

Here $Q_k$ is an orthonormal basis for $\mathcal K_k(A,q_1)$.
:::

:::{prf:proof}
The vector $v$ belongs to $\mathcal K_k$. Expanding its error gives

$$
v-x_1
=\frac1{c_1}X_2
\operatorname{diag}\bigl(p(\lambda_2),\ldots,p(\lambda_n)\bigr)c_{\mathrm{rest}}.
$$

Orthogonal projection gives the closest vector in the subspace, so

$$
\begin{aligned}
\|(I-Q_kQ_k^H)x_1\|_2
&\le\|x_1-v\|_2\\
&\le\frac{\|X_2\|_2\|c_{\mathrm{rest}}\|_2}{|c_1|}
\max_{i>1}|p(\lambda_i)|.
\end{aligned}
$$
:::

The first inequality expresses the best-approximation property of orthogonal projection. The second tells us how accurately the chosen polynomial isolates $x_1$. A small starting coefficient $c_1$ makes the bound worse. For a normal matrix, the eigenvector basis can be unitary and $\|X_2\|_2=1$.

For the Ritz-root filter, substitute $p(z)=s_\ell(z)/s_\ell(\lambda_1)$ into the second inequality:

$$
\left\|x_1-\frac{p(A)q_1}{c_1}\right\|_2
\le
\frac{\|X_2\|_2\|c_{\mathrm{rest}}\|_2}{|c_1|}
\max_{i>1}
\left|\frac{s_\ell(\lambda_i)}{s_\ell(\lambda_1)}\right|.
$$

This now bounds a particular Ritz direction. With $v_\ell=s_\ell(A)q_1$ and its unit Ritz vector $u_\ell$ defined above,

$$
\begin{aligned}
\frac{p(A)q_1}{c_1}
&=\frac{s_\ell(A)q_1}{c_1s_\ell(\lambda_1)}\\
&=\frac{\|v_\ell\|_2}{c_1s_\ell(\lambda_1)}\,u_\ell.
\end{aligned}
$$

Thus the filtered vector is a scalar multiple of the Ritz vector. If the other Ritz roots make $s_\ell$ small at unwanted eigenvalues relative to its value at $\lambda_1$, the bound is small. Normalizing then gives an accurate Ritz direction, up to a sign or phase, by the same normalization argument used for the power method.

A simpler trial polynomial can help compare the approximation available in a Krylov space with the power method. For example, suppose $\lambda_1=10$ and all unwanted eigenvalues lie in $[1,2]$. The polynomial

$$
p(z)=\left(\frac{z-1.5}{8.5}\right)^{k-1}
$$

satisfies $p(10)=1$ and $|p(\lambda_i)|\le(1/17)^{k-1}$ for $i>1$. This gives a faster-decaying subspace error bound than the power-method choice $(z/10)^{k-1}$, whose corresponding bound is $(1/5)^{k-1}$. Arnoldi makes these polynomial combinations available without explicitly choosing this polynomial.

This example bounds the accuracy available in the space; the Ritz-root filter identifies a particular computed direction. To certify an actual eigenvalue approximation, we return to the residual obtained directly from Arnoldi.

### What a Small Residual Guarantees

For a computed Ritz pair, the residual check from the beginning of the page gives a quantity we can evaluate without knowing any exact eigenvalues. Recall that $H_ky=\theta y$, $\|y\|_2=1$, and $u=Q_ky$. Then

$$
\boxed{
\|r\|_2
=\|f_k\|_2\,|(y)_k|
=|h_{k+1,k}|\,|(y)_k|
\le |h_{k+1,k}|.
}
$$

The coefficient $h_{k+1,k}=\|f_k\|_2$ measures how much of $Aq_k$ lies outside the current subspace. The other products $Aq_1,\ldots,Aq_{k-1}$ already lie inside it. Thus, if $h_{k+1,k}$ is small, the whole subspace is nearly invariant, and every normalized Ritz pair has a small residual. At exact breakdown, $f_k=0$ and these pairs are exact.

An individual pair can also have a small residual because its last coordinate $(y)_k$ is small. We therefore monitor the product $|h_{k+1,k}|\,|(y)_k|$, rather than requiring $h_{k+1,k}$ itself to decrease. The next theorem converts this Arnoldi residual into an eigenvalue error bound.

:::{prf:theorem} Residual Bound for Eigenvalue Error
:label: thm:non_hermitian_convergence

Suppose $A=X\Lambda X^{-1}$ and $\|u\|_2=1$. For $r=Au-\theta u$,

$$
\boxed{
\begin{aligned}
\min_i|\theta-\lambda_i|
&\le\kappa_2(X)\|r\|_2,\\
\kappa_2(X)&=\|X\|_2\|X^{-1}\|_2.
\end{aligned}
}
$$

If $A$ is normal, we can take $X$ unitary and obtain

$$
\min_i|\theta-\lambda_i|\le\|r\|_2.
$$
:::

:::{prf:proof}
If $\theta$ is an eigenvalue, the bound is immediate. Otherwise,

$$
u=X(\Lambda-\theta I)^{-1}X^{-1}r.
$$

Taking norms and using $\|u\|_2=1$ gives

$$
1\le
\frac{\kappa_2(X)\|r\|_2}{\min_i|\lambda_i-\theta|},
$$

which proves the result.
:::

This is the residual form of the Bauer–Fike bound; see [Saad, Section 3.2](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf). It guarantees proximity to some eigenvalue, without identifying which one. For a nonnormal matrix, $\kappa_2(X)$ can be large, so a small residual may give a much weaker eigenvalue error guarantee.

For a Ritz value, substituting the Arnoldi residual gives the explicit bound

$$
\boxed{
\min_i|\theta-\lambda_i|
\le\kappa_2(X)|h_{k+1,k}|\,|(y)_k|.
}
$$

For normal $A$, the factor $\kappa_2(X)$ can be taken as one.

Even when eigenvalue sensitivity makes the preceding bound weak, a small residual still gives a precise backward-error guarantee: the Ritz pair is exact for a nearby matrix. This statement does not require diagonalizability. To see it, set

$$
E=-ru^H.
$$

Then

$$
(A+E)u=Au-r=\theta u,
\qquad
\|E\|_2=\|r\|_2.
$$

Moreover, any perturbation making this pair exact must satisfy $Eu=-r$, hence $\|E\|_2\ge\|r\|_2$. The residual norm is therefore the smallest perturbation to $A$, in the matrix 2-norm, that makes the unit-vector pair exact. This statement remains useful even when eigenvalues are sensitive to perturbations.

For symmetric matrices, an orthonormal eigenvector basis removes the eigenvector condition factor from the bounds, and the projected matrix is tridiagonal. The next two sections exploit this structure through [Lanczos iteration](lanczos_iteration.md) and its [convergence analysis](lanczos_convergence.md).
