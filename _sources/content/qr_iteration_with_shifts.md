# Shifted QR Iteration

The [previous section](upper_hessenberg.md) makes each QR step inexpensive by reducing $A$ to upper Hessenberg form. A **shift** addresses a different problem: how many steps are needed before an eigenvalue, or a pair of eigenvalues, can be separated from the rest of the matrix. A poor shift can leave the iteration cycling. A good shift can make the separating subdiagonal entry nearly zero in a single step.

We focus on a real, generally nonsymmetric matrix. Its eigenvalues may be real or occur in complex conjugate pairs. The target is therefore **real Schur form**: an upper block triangular matrix with $1\times1$ blocks for real eigenvalues and $2\times2$ blocks for nonreal conjugate pairs.

## The Shifted QR Step

For $A\in\mathbb R^{n\times n}$, start from the Hessenberg reduction

$$
H=Q_H^TAQ_H,
$$

where $Q_H$ is the orthogonal matrix from the previous section. Throughout this page, $H$ denotes the current Hessenberg matrix. Its size is $m\times m$, initially with $m=n$, and its entries are $h_{ij}$.

For one step, choose a real shift $\mu$ and factor $H-\mu I=UR$, where $U$ is orthogonal and $R$ is upper triangular. The updated matrix, denoted by $H_{\mathrm{new}}$, is obtained by reversing the factors and adding the shift back:

$$
\boxed{
\begin{aligned}
H-\mu I&=UR,\\
H_{\mathrm{new}}&=RU+\mu I.
\end{aligned}}
$$

We replace $H$ by $H_{\mathrm{new}}$ before taking the next step. The shift changes the QR factors, but the update still preserves the eigenvalues. Since $R=U^T(H-\mu I)$,

$$
\begin{aligned}
H_{\mathrm{new}}
&=U^T(H-\mu I)U+\mu I\\
&=U^THU-\mu I+\mu I\\
&=U^THU.
\end{aligned}
$$

The Hessenberg structure is preserved by the adjacent-rotation implementation from the previous section, so a step costs $O(m^2)$ operations. A complex shift generally introduces complex arithmetic; below we explain how to handle conjugate shifts together while keeping the computation real.

## Why a Good Shift Helps

For intuition, first hold $\mu$ fixed. Shifted QR iteration on $H$ is then ordinary QR iteration on $H-\mu I$, followed by adding $\mu I$ back. If the eigenvalues of $H$ are $\lambda_j$, those of $H-\mu I$ are $\lambda_j-\mu$.

Recall that orthogonal iteration captures the directions associated with larger eigenvalue magnitudes first. If $\mu$ is close to one eigenvalue $\lambda_\star$, then $|\lambda_\star-\mu|$ is small. The other directions occupy the leading columns, leaving this eigenvalue at the bottom.

More precisely, suppose $H$ is diagonalizable, $\mu$ is not an eigenvalue, and a simple eigenvalue $\lambda_\star$ is uniquely closest to $\mu$. Under the starting-space assumption of {prf:ref}`thm:orth-it-nh-rigorous`, the leading $(m-1)$-dimensional subspace converges with the asymptotic factor

$$
\rho_\mu=
\frac{|\lambda_\star-\mu|}
{\min_{\lambda_j\ne\lambda_\star}|\lambda_j-\mu|}
<1.
$$

The smaller this ratio, the faster the last eigenvalue separates. Closeness to the target must be judged relative to its distance from the other eigenvalues.

In practice, the shift is recomputed after each step. As the trailing part of the matrix separates, it can give a better eigenvalue estimate, which makes the next shift more effective. This feedback explains the rapid local convergence often observed, but does not guarantee convergence from every starting matrix.

### A Shift That Cycles, and One That Deflates

Consider the nonsymmetric matrix

$$
H=\begin{pmatrix}0&2\\1&0\end{pmatrix},
\qquad \lambda=\pm\sqrt2.
$$

The simplest shift is the bottom-right entry, so $\mu=0$. A QR factorization with positive diagonal entries in $R$ is

$$
H=
\underbrace{\begin{pmatrix}0&1\\1&0\end{pmatrix}}_U
\underbrace{\begin{pmatrix}1&0\\0&2\end{pmatrix}}_R.
$$

Reversing the factors gives $\begin{pmatrix}0&1\\2&0\end{pmatrix}$. Its bottom-right entry is again zero, and the next step returns to $H$:

$$
\begin{pmatrix}0&2\\1&0\end{pmatrix}
\ \longrightarrow\
\begin{pmatrix}0&1\\2&0\end{pmatrix}
\ \longrightarrow\
\begin{pmatrix}0&2\\1&0\end{pmatrix}.
$$

Neither eigenvalue separates. Both are equally far from the shift.

With the exact shift $\mu=\sqrt2$, a QR factorization instead gives

$$
H-\sqrt2 I
=
\underbrace{\frac1{\sqrt3}
\begin{pmatrix}-\sqrt2&1\\1&\sqrt2\end{pmatrix}}_U
\underbrace{\begin{pmatrix}\sqrt3&-\sqrt6\\0&0\end{pmatrix}}_R.
$$

Consequently,

$$
H_{\mathrm{new}}=RU+\sqrt2 I
=\begin{pmatrix}-\sqrt2&-1\\0&\sqrt2\end{pmatrix}.
$$

One step makes the matrix triangular. An approximate shift also helps: with $\mu=1.414$, the new subdiagonal entry is approximately $2.01\times10^{-4}$, compared with its original value $1$. The example shows both stagnation and rapid progress. For a larger matrix, we need an inexpensive way to obtain such shift estimates without already knowing its eigenvalues.

## Deflation: Recognizing That an Eigenvalue Has Separated

A successful shift makes a subdiagonal entry small enough to replace by zero. Because $H$ is Hessenberg, this separates it into an upper block triangular matrix. Its eigenvalues are the eigenvalues of its diagonal blocks, so the smaller problems can be solved separately. This reduction is called **deflation**.

For a single real eigenvalue, the separating entry is $h_{m,m-1}$. Setting it to zero gives

$$
H\approx
\begin{pmatrix}
H_{11}&*\\
0&h_{mm}
\end{pmatrix}.
$$

Here $H_{11}$ consists of the first $m-1$ rows and columns. The asterisk denotes the remaining entries above $h_{mm}$, which do not affect the eigenvalues of this block triangular matrix. The scalar $h_{mm}$ is an eigenvalue of the matrix after this small perturbation. We record it and continue QR iteration on $H_{11}$. This remaining matrix is called the **active block**: it contains the eigenvalues still to be computed. We keep the notation $H$ for it and decrease $m$ by one.

For a trailing conjugate pair, the separating entry is instead $h_{m-1,m-2}$, when $m\ge3$. Setting it to zero isolates the trailing $2\times2$ block $B$:

$$
H\approx
\begin{pmatrix}
H_{11}&*\\
0&B
\end{pmatrix}.
$$

Now $H_{11}$ consists of the first $m-2$ rows and columns. The internal entry $h_{m,m-1}$ of $B$ need not become small. Once $B$ is isolated, an orthogonal $2\times2$ similarity puts it in real Schur form: triangular for real roots, or a $2\times2$ block for a nonreal pair. We compute its two eigenvalues directly and continue on $H_{11}$, again denoting this smaller matrix by $H$, with $m$ decreased by two. An isolated nonreal block is a completed part of the calculation, not a convergence failure. The same splitting principle applies to small subdiagonal entries anywhere inside the active block.

How small is small enough? Setting $h_{j+1,j}$ to zero changes $H$ by a matrix with 2-norm $|h_{j+1,j}|$. This gives a backward-error interpretation of the deflation test. A common initial check compares it with nearby diagonal entries:

$$
|h_{j+1,j}|\le
\tau\bigl(|h_{jj}|+|h_{j+1,j+1}|\bigr),
$$

where $\tau$ is of the order of machine precision. This is only an initial check. Robust nonsymmetric implementations use additional tests involving the off-diagonal entries, together with underflow safeguards. A small discarded entry gives a small matrix perturbation; for a sensitive eigenvalue problem, it need not give an equally small relative eigenvalue error.

If Schur vectors or the full Schur form are required, the active-block transformations must also update the rest of the matrix and the accumulated orthogonal factor. Deflation reduces the part still needing iteration; it does not mean discarding information needed for those outputs.

## Computing Shifts from the Trailing Block

All shift choices below use the current matrix $H$, whose size decreases as eigenvalues are deflated.

### The Bottom-Right Entry

The choice

$$
\mu=h_{mm}
$$

is often called the **Rayleigh quotient shift**, since $h_{mm}=\boldsymbol e_m^TH\boldsymbol e_m$, where $\boldsymbol e_m$ is the last coordinate vector. When the last row is nearly separated, this entry approximates a real eigenvalue. It costs nothing to compute, but the cycling example shows why it cannot be the only shift rule.

The rapid local convergence can be seen in a $2\times2$ calculation. Write the entry below the diagonal as $\varepsilon$ and choose the bottom-right entry as the shift:

$$
H=\begin{pmatrix}a&b\\\varepsilon&d\end{pmatrix},
\qquad \mu=d.
$$

Let $r=\sqrt{(a-d)^2+\varepsilon^2}$ be the length of the first column of $H-dI$, and assume $r>0$. A rotation-based QR factorization is

$$
H-dI=
\underbrace{\frac1r
\begin{pmatrix}a-d&-\varepsilon\\\varepsilon&a-d\end{pmatrix}}_U
\underbrace{\begin{pmatrix}
r&b(a-d)/r\\
0&-b\varepsilon/r
\end{pmatrix}}_R.
$$

The updated matrix is $H_{\mathrm{new}}=RU+dI$. Denote its entry below the diagonal by $\varepsilon_{\mathrm{new}}$. Multiplying the last row of $R$ by the first column of $U$ gives

$$
\begin{aligned}
\varepsilon_{\mathrm{new}}
&=\left(-\frac{b\varepsilon}{r}\right)
  \frac{\varepsilon}{r}\\
&=-\frac{b\varepsilon^2}{(a-d)^2+\varepsilon^2}.
\end{aligned}
$$

When $|\varepsilon|\ll|a-d|$, the denominator is approximately $(a-d)^2$, so

$$
|\varepsilon_{\mathrm{new}}|
\approx\frac{|b|}{(a-d)^2}|\varepsilon|^2.
$$

This is **quadratic local convergence**: the new coupling is proportional to the square of the old one. With $a,b,d$ fixed and $a\ne d$, a sufficiently small coupling becomes much smaller after one step. This calculation does not explain how to escape the cycle above, where $a=d=0$.

### The Trailing $2\times2$ Block

A better estimate uses the last two coordinates together:

$$
B=\begin{pmatrix}a&b\\c&d\end{pmatrix}
=\begin{pmatrix}
h_{m-1,m-1}&h_{m-1,m}\\
h_{m,m-1}&h_{mm}
\end{pmatrix}.
$$

Its two eigenvalues solve

$$
z^2-sz+t=0,
\qquad s=a+d,\qquad t=ad-bc.
$$

Equivalently, with

$$
\alpha=\frac{a+d}{2},
\qquad
\Delta=\left(\frac{a-d}{2}\right)^2+bc,
$$

the roots are

$$
\mu_1,\mu_2=
\begin{cases}
\alpha\pm\sqrt\Delta,&\Delta\ge0,\\
\alpha\pm i\sqrt{-\Delta},&\Delta<0.
\end{cases}
$$

When the roots are real, one possible single shift is the root closest to $d=h_{mm}$, with a fixed rule to break ties. This uses information that the bottom-right entry alone misses. In the cycling example, it selects one of $\pm\sqrt2$ and ends the cycle immediately.

When the roots are nonreal, they estimate a conjugate pair. A real matrix cannot become upper triangular over the reals while retaining nonreal eigenvalues. We must separate a two-dimensional block instead of trying to force its two eigenvalues onto real diagonal entries.

The roots of $B$ are estimates of eigenvalues of $H$; they need not be accurate until the trailing block is nearly separated. The displayed formulas explain the choices, but a numerical implementation should use scaled formulas that avoid overflow and cancellation. LAPACK's [DLANV2 routine](https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dlanv2.f) computes a real $2\times2$ Schur decomposition and its eigenvalues with such safeguards.

## Two Shifts in Real Arithmetic: The Francis Step

Applying a complex shift and then its conjugate would return, up to QR phase conventions, to a real problem. The **Francis double-shift step** combines the pair without carrying out the intervening complex computation.

For shifts $\mu_1$ and $\mu_2$ that are real or complex conjugates, form the polynomial

$$
p(z)=(z-\mu_1)(z-\mu_2)=z^2-sz+t,
\qquad
s=\mu_1+\mu_2,\quad t=\mu_1\mu_2.
$$

Both $s$ and $t$ are real. For the eigenvalues of $B$, these coefficients are simply $s=\operatorname{tr}(B)$ and $t=\det(B)$. Thus

$$
p(H)=H^2-sH+tI
$$

is real even when neither shift is real. A conceptual double-shift step takes the QR factorization of $p(H)$. Using $U$ and $R$ for these factors, the update is

$$
\boxed{
p(H)=UR,\qquad H_{\mathrm{new}}=U^THU.
}
$$

Here $U$ is orthogonal, $R$ is upper triangular, and $H_{\mathrm{new}}$ is the matrix after the double-shift step. When $p(H)$ is nonsingular, this combines two successive shifted QR steps, up to the signs or phases allowed in their QR factors.

To see why, first apply the shift $\mu_1$. Denote the intermediate matrix by $H_1$:

$$
H-\mu_1I=U_1R_1,\qquad
H_1=R_1U_1+\mu_1I=U_1^H H U_1.
$$

These factors may be complex, so the superscript $H$ denotes the conjugate transpose. Apply the second shift $\mu_2$ to $H_1$, obtaining $H_2$:

$$
H_1-\mu_2I=U_2R_2,\qquad
H_2=R_2U_2+\mu_2I=U_2^H H_1U_2.
$$

The combined change of basis is $U_1U_2$, since

$$
H_2=(U_1U_2)^H H(U_1U_2).
$$

This same product appears in a QR factorization of $p(H)$:

$$
\begin{aligned}
p(H)
&=(H-\mu_2I)(H-\mu_1I)\\
&=(H-\mu_2I)U_1R_1\\
&=U_1(H_1-\mu_2I)R_1\\
&=(U_1U_2)(R_2R_1).
\end{aligned}
$$

The product $R_2R_1$ is upper triangular, so $U_1U_2$ is a unitary QR factor of the real matrix $p(H)$. The real factor $U$ in $p(H)=UR$ differs from $U_1U_2$ only by column signs or unit complex factors. Thus $H_{\mathrm{new}}$ carries out the same two-shift transformation, up to those harmless choices of basis.

The efficient implementation of this step requires the implicit QR property, which is beyond the scope of this section.

The trailing-block pair is a standard choice, but implementations can modify it. For example, when the two roots are real, LAPACK's [DLAHQR routine](https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dlahqr.f) uses the root closer to $d$ twice. Conjugate roots are kept together. Double shifting handles complex pairs in real arithmetic; it does not by itself rule out stagnation.

## Stagnation and Practical Remedies

The trailing block is useful because it is small and becomes informative near convergence. Far from convergence, it may keep producing ineffective shifts. A robust implementation therefore monitors progress and changes strategy when needed.

### Exceptional Shifts

If many sweeps pass without deflation, repeating the same trailing-block rule can sustain a cycle or a nearly stationary pattern. Here a **sweep** means one pass of the QR transformations through the active block. **Exceptional shifts** deliberately depart from the usual shift rule to break the pattern. They are usually constructed from local matrix entries and subdiagonal magnitudes, so they adapt to the matrix's scale. They need not be better eigenvalue estimates; their purpose is to move the iteration away from stagnation.

For example, LAPACK's DLAHQR counts iterations since the last deflation and periodically substitutes exceptional shifts. It also uses conservative deflation tests and an iteration limit. These safeguards are part of the method, rather than an assumption that every trailing-block shift will work.

### Larger Trailing Windows

Closely spaced eigenvalues may be difficult to separate, and for a strongly nonnormal matrix the trailing $2\times2$ eigenvalues can be poor estimates of the full spectrum. A **trailing window** is a block formed from the last several rows and columns of $H$. Using a larger window can supply several shift estimates. **Multishift QR iteration** applies these shifts in a structured sweep.

Modern implementations also use **aggressive early deflation**: compute a small Schur decomposition of a trailing window and test whether some of its Schur blocks have negligible coupling to the rest of the matrix. This can recognize converged eigenvalues before the ordinary subdiagonal tests do. LAPACK's [implementation notes for xHSEQR](https://www.netlib.org/lapack/lawnspdf/lawn187.pdf) describe how larger windows, multiple shifts, and early deflation are combined, including a fallback when double-shift iteration fails to converge.

These strategies improve reliability, but the usual trailing-block shift rules do not provide an unconditional convergence guarantee for every nonsymmetric matrix. A solver must impose an iteration limit and report any unresolved block. LAPACK's [DHSEQR documentation](https://www.netlib.org/lapack/double/dhseqr.f) specifies both successful real Schur output and partial results on nonconvergence.
