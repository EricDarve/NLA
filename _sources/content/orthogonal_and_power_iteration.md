# Power and Orthogonal Iteration Methods

Given the way the orthogonal iteration method was constructed, we should expect some connection to the power iteration method. We explore this in more detail here. We show how the sequence of matrices $Q_k$ from the orthogonal iteration method corresponds to a QR factorization of the sequence of matrices $A^k$ from the power iteration method.

(sec:oi-qr-ak)=
## $Q_k$ as the $Q$-factor in a QR factorization of $A^k$

Let $A\in\mathbb{C}^{n\times n}$. Orthogonal iteration with an initial matrix $Q_0\in\mathbb{C}^{n\times p}$ whose columns are orthonormal proceeds by computing, at each step, the thin QR factorization

$$
AQ_k = Q_{k+1} R_{k+1}
$$

where $Q_k$ is unitary and $R_k$ is upper triangular.

Define the accumulated upper triangular factor

$$
\widehat{R}_k \;:=\; R_k R_{k-1}\cdots R_1 \quad (k\ge 1),\qquad \widehat{R}_0 := I_p.
$$

:::{prf:theorem} QR of $A^k Q_0$
:label: thm:AkQR

For all $k\ge 0$,

$$
A^k\,Q_0 \;=\; Q_k\,\widehat{R}_k .
$$

In particular, if $Q_0=I_n$ (full orthogonal iteration on all of $\mathbb{C}^n$), then

$$
A^k \;=\; Q_k\,\widehat{R}_k ,
$$

so $Q_k$ is the $Q$-factor and $\widehat{R}_k$ the $R$-factor in the QR factorization of $A^k$.
:::

:::{prf:proof}
*Base cases.* 

For $k=0$, the identity is $Q_0 = Q_0$. 

For $k=1$, it reads $A Q_0 = Q_1 R_1$, which is the defining QR step.

*Inductive step.* Assume $A^k Q_0 = Q_k \widehat{R}_k$. Then

$$
\begin{aligned}
A^{k+1}Q_0
&= A(A^k Q_0)
= A(Q_k \widehat{R}_k)
= (A Q_k)\,\widehat{R}_k \\
&= (Q_{k+1} R_{k+1})\,\widehat{R}_k
= Q_{k+1}\,(R_{k+1}\widehat{R}_k)
= Q_{k+1}\,\widehat{R}_{k+1},
\end{aligned}
$$

which proves the claim by induction.
:::

:::{admonition} Remarks
:class: note

- A product of upper triangular matrices is upper triangular, and with positive diagonals in each $R_j$, the factor $\widehat{R}_k$ is upper triangular with positive diagonal. Hence the factorization above is the unique thin QR factorization of $A^k Q_0$.
- When $p=1$ (the power method), $\widehat{R}_k$ is the scalar $\prod_{j=1}^k r_j$ with $r_j=\lVert A q_{j-1}\rVert_2$, and $q_k$ is the normalized power-iteration vector. The identity reduces to the familiar normalization of $A^k q_0$.
- The column-space relation $\operatorname{range}(Q_k)=\operatorname{range}(A^k Q_0)$ shows that orthogonal iteration orthonormalizes the images of the initial subspace under repeated application of $A$. This is the subspace analogue of the power method and explains why the nested leading-column spans of $Q_k$ converge to invariant subspaces under standard spectral gap assumptions.
:::