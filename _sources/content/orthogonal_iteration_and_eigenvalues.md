# Orthogonal Iteration and Eigenvalues

The [previous section](orthogonal_iteration.md) establishes how orthogonal iteration finds invariant subspaces and Schur directions. To use these results for eigenvalue computation, we represent the action of $A$ in the current basis.

We keep the same notation: $Q_\star$ is an exact Schur basis in $A=Q_\star TQ_\star^H$, with $|\lambda_1|\ge\cdots\ge|\lambda_n|$, and $Q_k\in\mathbb C^{n\times p}$ is the orthonormal basis at iteration $k$.

## Extracting Eigenvalues from a Subspace

For $p<n$, a vector in the current subspace has the form $u=Q_ky$. Multiplication by $A$ can take it outside the subspace. Projecting back gives

$$
T_k=Q_k^HAQ_k,
\qquad
Q_kQ_k^HAu=Q_kT_ky.
$$

The $p\times p$ matrix $T_k$ therefore describes the projected action of $A$. For $p=1$, it is the Rayleigh quotient used in the power method. For general $p$, its eigenvalues are called **Ritz values**. Take a unit eigenvector of $T_k$:

$$
T_ky=\theta y,
\qquad \|y\|_2=1.
$$

The corresponding **Ritz vector** $u=Q_ky$ also has unit norm. Using the subspace residual $E_k=AQ_k-Q_kT_k$ from the previous section,

$$
\begin{aligned}
Au-\theta u
&=AQ_ky-Q_kT_ky\\
&=E_ky,\\
\|Au-\theta u\|_2
&\le\|E_k\|_2.
\end{aligned}
$$

Because $Q_k^HE_k=Q_k^HAQ_k-T_k=0$, this residual is perpendicular to the current subspace. Thus the Ritz pair satisfies the eigenvector equation after projection. If the subspace is exactly invariant, $E_k=0$ and $(\theta,u)$ is an exact eigenpair of $A$. For a nonnormal matrix, the eigenvalue error can be much larger than the residual.

## From Subspace Convergence to Schur Form

### One Gap: Convergence of the Subspace

For $p<n$, suppose the preceding convergence theorem applies: $A=X\Lambda X^{-1}$ is diagonalizable, $|\lambda_p|>|\lambda_{p+1}|$, and the first $p$ rows of $X^{-1}Q_0$ form an invertible matrix. The columns of $X$ follow the eigenvalue ordering above.

This gap controls separation of the wanted subspace from the unwanted directions. No gaps within the wanted group are assumed. Thus $\operatorname{range}(Q_k)$ approaches $\operatorname{range}(Q_{\star,p})$, but **the basis matrix $Q_k$ need not converge**. Two orthonormal vectors can keep rotating within a plane even though the plane itself stays fixed.

To compare the current basis with the exact one, we need a general unitary $p\times p$ matrix $W_k$ that can mix its columns. Subspace convergence guarantees that such matrices exist with

$$
\|Q_kW_k^H-Q_{\star,p}\|_2\longrightarrow0.
$$

The alignment depends on the unknown target basis $Q_{\star,p}$, so $W_k$ is not available as a step of the algorithm. Its role here is theoretical: it lets us express what subspace convergence implies for $T_k$. Writing the leading $p\times p$ Schur block as

$$
T_{\star,p}=Q_{\star,p}^HAQ_{\star,p},
$$

we obtain

$$
\|W_kT_kW_k^H-T_{\star,p}\|_2\longrightarrow0.
$$

Unitary changes of basis preserve eigenvalues, and eigenvalues depend continuously on matrix entries. Thus the Ritz values approach $\lambda_1,\ldots,\lambda_p$, up to reordering and counting multiplicities. However, this argument does not give a computable rotation that makes $T_k$ triangular, nor does it show that its diagonal entries approach the eigenvalues.

### Strict Ordering: Convergence of Schur Directions

For full iteration, $p=n$, $T_k$ is unitarily similar to $A$. This guarantees that the two matrices have the same eigenvalues in exact arithmetic, but it does not compute their values. Finding the eigenvalues of a general $T_k$ is still the original eigenvalue problem. Convergence of the individual Schur directions provides a way to read the eigenvalues from its diagonal.

As proved on the previous page, this stronger conclusion holds when

$$
|\lambda_1|>\cdots>|\lambda_n|>0,
$$

and the starting-space condition holds for every leading set of columns. Now every leading subspace converges, fixing each Schur direction. The only remaining freedom is a sign or unit complex phase for each column. Following the power method's notation $\phi_k$ for a single phase, write these phases as $\phi_{k,j}$ and collect them in a **diagonal** unitary matrix

$$
\Phi_k=\operatorname{diag}(\phi_{k,1},\ldots,\phi_{k,n}),
\qquad |\phi_{k,j}|=1.
$$

Whereas $W_k$ may mix columns, $\Phi_k$ only adjusts their phases, giving

$$
\begin{aligned}
\|Q_k\Phi_k^H-Q_\star\|_2&\longrightarrow0,\\
\|\Phi_kT_k\Phi_k^H-T\|_2&\longrightarrow0.
\end{aligned}
$$

We do not need to compute $\Phi_k$: these phase adjustments leave diagonal entries unchanged and preserve the magnitudes of all other entries. Since $T$ is upper triangular, the unadjusted matrix $T_k$ therefore satisfies

$$
\boxed{(T_k)_{jj}\longrightarrow\lambda_j,
\qquad (T_k)_{j\ell}\longrightarrow0\quad(j>\ell).}
$$

This gives a way to compute the eigenvalues from the diagonal of $T_k$ without knowing $Q_\star$ or correcting the phases. Entries above the diagonal may keep changing with those phases.

## From Orthogonal Iteration to QR Iteration

The diagonal approximations can also guide acceleration through **shifts**. Replacing $A$ by $A-\mu I$ preserves its invariant subspaces while changing the eigenvalue ratios that control convergence. Suitable choices can speed up convergence; their selection and analysis are developed in [shifted QR iteration](qr_iteration_with_shifts.md).

To turn these results into a practical eigenvalue method, the [next section](qr_iteration.md) derives an update from $T_k$ directly to $T_{k+1}$ using a QR factorization. This **QR iteration** produces the same matrix sequence without explicitly updating the full basis $Q_k$ and forming $Q_k^HAQ_k$ at every step.
