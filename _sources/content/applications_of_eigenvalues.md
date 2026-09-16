# Applications of Eigenvalues

Eigenvectors identify directions or patterns that a matrix acts on by scaling. The interpretation of the eigenvalues depends on the model: they can describe growth rates, squared vibration frequencies, variances, or measurable energies. The examples below show how these quantities arise.

## Connections to the Previous Sections

We have already seen how eigenvalues describe [time evolution](eigendecomposition.md): a mode evolves through $\lambda^k$ in discrete time and $e^{\lambda t}$ in continuous time. The heat-exchange example illustrates how this separates a preserved quantity from a decaying difference.

We also established, with eigenvalues counted with algebraic multiplicity,

$$\det(A)=\prod_{j=1}^n\lambda_j,\qquad
\operatorname{tr}(A)=\sum_{j=1}^n\lambda_j.$$

These identities hold even for defective matrices; their proof uses [Schur form](eigendecomposition.md). They connect eigenvalues to familiar matrix properties, without requiring us to compute eigenvalues just to find a determinant or trace.

## Vibrations: Mode Shapes and Natural Frequencies

Consider small displacements $u(t)$ of a system of unit masses connected by springs. An undamped linear model has the form

$$u''(t)+Ku(t)=0,$$

where $K$ is the **stiffness matrix**, describing the restoring forces. Assume $K$ is real symmetric positive definite. Look for a motion with a fixed spatial pattern $q$:

$$u(t)=q\cos(\omega t).$$

Substitution gives

$$Kq=\omega^2q.$$

Thus the eigenvector $q$ describes a **mode shape**, and its eigenvalue is the **square of the natural angular frequency**. If $\lambda>0$, then $\omega=\sqrt{\lambda}$, measured in radians per unit time.

For two unit masses connected to each other and to fixed walls by three springs of unit stiffness,

$$K=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}.$$

The mode $(1,1)^T$ has eigenvalue $1$: the masses move together with angular frequency $1$. The mode $(1,-1)^T$ has eigenvalue $3$: they move in opposite directions with angular frequency $\sqrt3$.

For a general mass matrix $M$, the equation is $Mu''+Ku=0$, and the mode equation becomes

$$Kq=\omega^2Mq.$$

This is a **generalized eigenvalue problem**. External forcing near a natural frequency can produce large oscillations if it excites the corresponding mode. This is resonance; its response also depends on damping and the applied force.

## Principal Component Analysis: Directions of Greatest Variation

Suppose $z_1,\ldots,z_N\in\mathbb{R}^d$ are data vectors after subtracting their mean. Define the covariance matrix using the averaging convention

$$C=\frac1N\sum_{i=1}^N z_i z_i^T.$$

For a unit vector $q$, the scalar $q^Tz_i$ is the coordinate of observation $i$ along $q$. Its variance is

$$\frac1N\sum_{i=1}^N(q^Tz_i)^2=q^TCq.$$

The matrix $C$ is symmetric positive semidefinite. Write its orthonormal eigendecomposition as $C=Q\Lambda Q^T$, with $\lambda_1\ge\cdots\ge\lambda_d\ge0$. For $q=Qc$ with $\|c\|_2=1$,

$$q^TCq=\sum_{j=1}^d\lambda_jc_j^2\le\lambda_1.$$

Equality is attained by a leading eigenvector. This is the first **principal component direction**: the direction along which the data vary most. Successive eigenvectors maximize the remaining variance subject to orthogonality to the earlier directions.

Keeping the first $r$ eigenvectors gives the approximation

$$\widehat z_i=Q_rQ_r^Tz_i.$$

It preserves variance $\sum_{j=1}^r\lambda_j$, while the average squared reconstruction error is

$$\frac1N\sum_{i=1}^N\|z_i-\widehat z_i\|_2^2
=\sum_{j=r+1}^d\lambda_j.$$

For example, if $C=\begin{pmatrix}5&4\\4&5\end{pmatrix}$, its eigenvalues are $9$ and $1$. Projecting onto the leading direction $(1,1)^T/\sqrt2$ retains $90\%$ of the total variance. This quantifies what is retained; large variance need not mean importance for every scientific question. The [SVD](singular_value_decomposition.md) provides another way to describe this approximation.

## Quantum Mechanics: Measurement Outcomes

In a finite-dimensional quantum model, a measurable quantity is represented by a Hermitian matrix. For energy, this matrix is the **Hamiltonian** $H$. Its eigenvalue equation is

$$Hq_j=E_jq_j.$$

The real numbers $E_j$ are the possible measured energies, and the eigenvectors are corresponding **energy eigenstates**.

A normalized state $\psi$ can be expanded in an orthonormal eigenvector basis:

$$\psi=\sum_j c_jq_j,\qquad \sum_j|c_j|^2=1.$$

For a nonrepeated eigenvalue $E_j$, the probability of measuring $E_j$ is $|c_j|^2$. If an energy is repeated, add the squared coefficients in its eigenspace. Equivalently,

$$\Pr(E)=\|P_E\psi\|_2^2,$$

where $P_E$ is the orthogonal projection onto the eigenspace for energy $E$.

For example, let $H=\operatorname{diag}(0,2)$ and $\psi=(\sqrt3/2,1/2)^T$. The measured energy is $0$ with probability $3/4$ and $2$ with probability $1/4$. Here eigenvalues describe possible outcomes, while projections onto eigenspaces determine their probabilities.

## PageRank: A Stationary Distribution on a Network

PageRank models a user moving among web pages. Let $L_{ij}$ be the probability of following a link from page $j$ to page $i$. Divide each page's outgoing links equally among its destinations. If a page has no outgoing links, assign a uniform distribution over all $n$ pages. Thus every entry is nonnegative and every column sums to one.

To allow movement between otherwise disconnected parts of the network, suppose the user follows this link rule with probability $\alpha$, where $0<\alpha<1$, and otherwise jumps to a uniformly chosen page. The resulting transition matrix is

$$G=\alpha L+\frac{1-\alpha}{n}\mathbf{1}\mathbf{1}^T,$$

where $\mathbf{1}$ is the vector of all ones. For a column vector $p_k$ of page probabilities, $p_{k+1}=Gp_k$.

The PageRank vector is the stationary probability distribution:

$$Gp=p,\qquad p_i\ge0,\qquad \sum_i p_i=1.$$

Thus it is an eigenvector for eigenvalue $1$. The positive jump probability makes this stationary distribution unique and ensures convergence to it from any initial probability distribution. Its entries give the long-run fractions of visits to the pages, which serve as ranking scores.

Numerical methods for finding such eigenvectors will be discussed later in the course.
