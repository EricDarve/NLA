# CG vs. GMRES

The [conjugate gradient method (CG)](conjugate_gradient.md) and [GMRES](gmres.md) both approximate the solution of $Ax=b$ in growing Krylov spaces. Their main differences are the quantity they minimize and the work needed to build each new direction.

Their basis constructions explain the difference in cost:

* **CG** is connected to the **Lanczos process**, which uses symmetry to obtain a short recurrence.
* **GMRES** uses the **Arnoldi process**, which constructs an orthonormal basis for a general square matrix.

Symmetry gives the short recurrence; positive definiteness gives CG its energy norm and minimization property.

## Key Differences

### Applicability

* **CG:** The standard theory assumes that $A$ is **symmetric positive definite (SPD)**. For an indefinite matrix, the step denominator $p_k^TAp_k$ can vanish, so the same guarantees do not apply.
* **GMRES:** It applies to **any square, nonsingular matrix** $A$, including nonsymmetric and symmetric indefinite matrices. It also works for SPD matrices, although CG usually requires less work and storage.

### Recurrence and Optimality

For an initial guess $x_0$, let $r_0=b-Ax_0$ and $\mathcal K_k=\mathcal K_k(A,r_0)$. Both methods search in $x_0+\mathcal K_k$, but use different criteria:

$$
\begin{aligned}
\text{CG:}\qquad
x_k&=\underset{z\in x_0+\mathcal K_k}{\operatorname{argmin}}\|x-z\|_A,\\
\text{GMRES:}\qquad
x_k&=\underset{z\in x_0+\mathcal K_k}{\operatorname{argmin}}\|b-Az\|_2.
\end{aligned}
$$

Here $x=A^{-1}b$ and $\|v\|_A=\sqrt{v^TAv}$. As shown in the CG derivation, its minimization can be performed without knowing $x$.

* **CG** uses short updates for the solution, residual, and search direction. A-orthogonality makes projections onto older search directions vanish, so only a few vectors are needed.
* **Standard GMRES** uses the full Arnoldi basis. For a general nonsymmetric matrix, each new vector must be orthogonalized against all preceding basis vectors $q_1,\ldots,q_k$.

In exact arithmetic, CG's energy error decreases, while GMRES's residual norm cannot increase. CG's Euclidean residual norm need not decrease at every step. These are different guarantees, so comparing residual histories alone does not compare their respective optimality properties.

### Computational and Storage Cost

For a sparse $n\times n$ matrix, let $\operatorname{nnz}(A)$ denote its number of nonzero entries. The costs below exclude storing $A$ and assume no preconditioning. CG's work per step stays fixed, while the orthogonalization work of unrestarted GMRES grows with the basis.

* **Conjugate Gradients (CG):**
    * **Cost per Iteration $k$:** $O(\operatorname{nnz}(A) + n)$, independent of $k$. It consists of one sparse matrix-vector product and a few vector operations.
    * **Storage:** $O(n)$. CG only needs to store a few vectors (the solution, residual, and direction vectors).

* **GMRES:**
    * **Cost per Iteration $k$:** $O(\operatorname{nnz}(A) + kn)$. The $kn$ term comes from orthogonalizing against the $k$ previous basis vectors and can dominate when $k$ is large.
    * **Storage:** $O(kn)$. GMRES stores the basis $Q_k=[q_1,\ldots,q_k]$ and the next vector. The small Hessenberg and QR factors require $O(k^2)$ additional storage, included in $O(kn)$ since $k\le n$.

### Comparison Summary Table

| Algorithm | Applicability | Work at Iteration $k$ | Working Storage | Total Work after $k$ Steps |
| :--- | :--- | :---: | :---: | :---: |
| **CG** | SPD | $O(\operatorname{nnz}(A) + n)$ | $O(n)$ | $O(k(\operatorname{nnz}(A) + n))$ |
| **GMRES** | Square, nonsingular | $O(\operatorname{nnz}(A) + kn)$ | $O(kn)$ | $O(k\operatorname{nnz}(A) + k^2n)$ |

### The Practicality of Restarting

The growing basis can make long GMRES runs expensive in both work and memory. **Restarted GMRES**, denoted GMRES($m$), limits this growth by keeping at most $m$ search directions per cycle.

1.  We choose a fixed restart parameter $m$ (e.g., $m=20$ or $m=50$).
2.  We run up to $m$ GMRES steps, stopping earlier if the requested accuracy is reached.
3.  We form the current approximation and use it as the new starting guess.
4.  We discard the basis and build a new one from the current residual.

Restarting limits **working storage to $O(mn)$** and **work per cycle to $O(m\operatorname{nnz}(A)+m^2n)$**. Total work also depends on how many cycles are needed.

The current approximation is retained, but the earlier search directions are discarded. Thus restarted GMRES no longer minimizes the residual over the full Krylov space accumulated since the original start. Restarting can require more iterations and can cause stagnation, as the example at the end of the [GMRES section](gmres.md) shows. CG's short recurrence already keeps storage fixed, so it needs no restart to control basis storage.

For symmetric matrices, we can retain the GMRES residual objective while avoiding its growing basis storage. The [next section](minres.md) develops this specialization, MINRES, which also applies to symmetric indefinite systems.
