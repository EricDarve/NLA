# A Comparative Summary: CG vs. GMRES

Both the Conjugate Gradients (CG) and GMRES algorithms are foundational Krylov subspace methods for solving large, sparse linear systems $Ax=b$. However, they are designed for different classes of problems, and their computational costs and practical behaviors differ dramatically.

The primary difference stems from the underlying process they use to build their respective Krylov subspaces:

* **CG** is built upon the **Lanczos process**, which only applies to symmetric matrices.
* **GMRES** is built upon the **Arnoldi process**, which applies to any square matrix.

This single distinction leads to all the major differences in their performance, cost, and applicability.

## Key Differences

### Applicability

* **Conjugate Gradients (CG):** This is a specialized algorithm. It is guaranteed to converge *only* if the matrix $A$ is **symmetric positive definite (SPD)**. It leverages this symmetry to achieve high efficiency.
* **GMRES:** This is a general-purpose algorithm. It can be applied to **any square, invertible matrix** $A$, including those that are non-symmetric, indefinite, or both.

### Recurrence and Optimality

* **Conjugate Gradients (CG):** By leveraging the $A$-norm, CG satisfies an optimality condition using a **short-term (three-term) recurrence**. This means that to compute the next vector, it only needs information from the previous two steps.
* **GMRES:** Lacking symmetry, GMRES must use a **long-term recurrence** (the full Arnoldi process). To enforce orthogonality, each new basis vector $q_{k+1}$ must be explicitly orthogonalized against *all* preceding vectors $q_1, \dots, q_k$. It finds the solution in the subspace that minimizes the 2-norm of the residual.

### Computational and Storage Cost

This is the most significant practical difference. The short-term recurrence of CG makes it extremely cheap, while the long-term recurrence of GMRES makes its cost grow with each iteration.

* **Conjugate Gradients (CG):**
    * **Cost per Iteration $k$:** $O(\text{nnz}(A) + n)$. This cost is **constant**. It is dominated by one sparse matrix-vector product ($O(\text{nnz})$) and a few vector operations ($O(n)$).
    * **Storage:** $O(n)$. CG only needs to store a few vectors (the solution, residual, and direction vectors).

* **GMRES:**
    * **Cost per Iteration $k$:** $O(\text{nnz}(A) + kn)$. This cost **grows linearly with $k$**. The $\text{nnz}(A)$ term is the matrix-vector product, but the $kn$ term comes from orthogonalizing the new vector against the $k$ previous basis vectors. This $kn$ term quickly becomes the dominant cost.
    * **Storage:** $O(kn)$. GMRES must store the entire basis $Q_k = [q_1, \dots, q_k]$ to perform the orthogonalization, so its memory footprint grows at each step.

### Comparison Summary Table

| Algorithm | Applicability | Cost at Iteration $k$ | Total Storage | Total Time (after $k$ steps) |
| :--- | :--- | :---: | :---: | :---: |
| **CG** | Symmetric Positive Definite | $O(\text{nnz} + n)$ | $O(n)$ | $O(k \, (\text{nnz} + n))$ |
| **GMRES** | Any Square Matrix | $O(\text{nnz} + kn)$ | $O(kn)$ | $O(k \cdot \text{nnz} + k^2n)$ |

### The Practicality of Restarting

The comparison table makes the primary weakness of GMRES clear: **it struggles with large $k$**.

The $k^2n$ time complexity and $kn$ storage cost mean that running GMRES for many iterations is computationally prohibitive. It becomes unacceptably slow, and its memory requirement will eventually exceed the system's capacity.

This leads to a common practical compromise: **restarted GMRES**, often denoted as **GMRES($m$)**.

1.  We choose a fixed restart parameter $m$ (e.g., $m=20$ or $m=50$).
2.  We run $m$ steps of the "inner" GMRES algorithm.
3.  We compute the approximate solution $x^{(m)}$ and use it as the new starting guess.
4.  We "restart" the algorithm, discarding the entire Krylov basis $Q_m$ and building a new one from the new residual $r^{(m)} = b - A x^{(m)}$.

This "restarting" procedure effectively caps the computational and storage costs at $O(mn)$ and $O(m \cdot \text{nnz} + m^2n)$, respectively.

The drawback, however, is significant: **restarting slows down convergence**. By throwing away the Krylov subspace, we lose all the information the algorithm has gathered about the matrix. This can lead to a much higher total number of iterations and, in some difficult cases, can cause the algorithm to stagnate completely. CG, by contrast, never needs to be restarted.