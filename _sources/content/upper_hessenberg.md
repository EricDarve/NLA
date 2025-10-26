# Reduction to Hessenberg Form

The standard QR iteration algorithm, as we have discussed, is a powerful theoretical tool. However, in its direct form, it is computationally prohibitive for large matrices. A naive QR factorization of a dense $n \times n$ matrix $A$ requires $O(n^3)$ floating-point operations (flops). If the iteration takes $k$ steps to converge, the total cost would be $O(k \cdot n^3)$, which is not competitive.

The key to making the QR iteration practical is to first reduce the matrix $A$ to a form that is "cheaper" to iterate on, without changing its eigenvalues. This is a classic two-phase approach:

1.  **Phase 1 (Reduction):** Transform $A$ into an **upper Hessenberg matrix** $H$ via a similarity transformation. This is a one-time cost of $O(n^3)$.
2.  **Phase 2 (Iteration):** Apply the QR iteration directly to $H$. Since $H$ has a special structure, this phase is much faster, costing only $O(n^2)$ per iteration.

The total cost becomes $O(n^3 + k \cdot n^2)$, a dramatic improvement over $O(k \cdot n^3)$.

### Why Upper Hessenberg Form?

An $n \times n$ matrix $H$ is **upper Hessenberg** if all its entries below the first subdiagonal are zero. That is, $h_{ij} = 0$ for all $i > j+1$.

$$
H = \begin{pmatrix}
h_{11} & h_{12} & h_{13} & \cdots & h_{1n} \\
h_{21} & h_{22} & h_{23} & \cdots & h_{2n} \\
0 & h_{32} & h_{33} & \cdots & h_{3n} \\
\vdots & \ddots & \ddots & \ddots & \vdots \\
0 & \cdots & 0 & h_{n,n-1} & h_{nn}
\end{pmatrix}
$$

This form is the ideal compromise for the QR iteration for two critical reasons:

1.  **It is computationally cheap to obtain.** As we will see, we can compute an orthogonal $Q$ such that $H = Q A Q^T$ in $O(n^3)$ flops. This is far cheaper than the iterative process required for a full **Schur decomposition** (which produces an upper triangular matrix).
2.  **It is invariant under the QR iteration.** If $H_k$ is upper Hessenberg, we compute its QR factorization $H_k = Q_k R_k$. The subsequent "reverse multiplication" step, $H_{k+1} = R_k Q_k$, also produces an upper Hessenberg matrix. This "preservation of form" is essential.
3.  **The QR factorization is fast.** Computing the QR factorization of an $n \times n$ upper Hessenberg matrix does not require $O(n^3)$ flops. It can be done in $O(n^2)$ time, for example, by using a sequence of $n-1$ **Givens rotations** to eliminate the subdiagonal elements.

### How the Reduction Works

Our goal is to find an orthogonal matrix $Q$ such that $H = Q A Q^T$ is upper Hessenberg. This is a **unitary similarity transformation**, which is critical because it **preserves the eigenvalues** of $A$.

We cannot introduce all the necessary zeros at once. Instead, we build $Q$ as a product of $n-2$ **Householder reflectors**, $Q = Q_{n-2} \cdots Q_2 Q_1$. The process works column by column.

**Step 1: The First Column**

We want to introduce zeros in the first column below the subdiagonal, i.e., in entries $(3, 1), (4, 1), \dots, (n, 1)$.

Find a Householder reflector, $Q_1$, that operates on rows 2 through $n$. This $Q_1$ is constructed to zero out the desired elements in the first column, exactly as in a standard Householder QR factorization.

$$
Q_1 A = \begin{pmatrix}
a_{11}' & a_{12}' & \cdots & a_{1n}' \\
a_{21}' & a_{22}' & \cdots & a_{2n}' \\
0 & a_{32}' & \cdots & a_{3n}' \\
\vdots & \vdots & \ddots & \vdots \\
0 & a_{n2}' & \cdots & a_{nn}'
\end{pmatrix}
$$

1.  Now, to maintain a similarity transformation, we *must* apply the reflector from the right: $A' = (Q_1 A) Q_1^T$. (Note: Householder reflectors are symmetric, so $Q_1^T = Q_1$).
2.  Because $Q_1$ only operates on rows/columns 2 through $n$, this right-multiplication $A' Q_1^T$ only mixes columns 2 through $n$. It **does not affect the first column**, so the zeros we just created are preserved.

**Step 2: The Second Column**

Our new matrix $A' = Q_1 A Q_1^T$ has the desired first column. We now move to the second column and aim to zero out entries $(4, 2), \dots, (n, 2)$.

1.  We find a new Householder reflector, $Q_2$, that operates on rows 3 through $n$.
2.  We apply it from the left and right: $A'' = Q_2 A' Q_2^T$.
3.  The left-multiplication $Q_2 A'$ creates the zeros in the second column.
4.  The right-multiplication $A'' Q_2^T$ only affects columns 3 through $n$. It does not destroy the zeros in the first column or the newly created zeros in the second column.

**Continuing the Process**

We repeat this process for a total of $n-2$ steps. At step $k$, we use a reflector $Q_k$ that operates on rows/columns $k+1$ through $n$ to introduce zeros in the $k$-th column.

The final matrix is:

$$
H = Q_{n-2} \cdots Q_2 Q_1 A Q_1^T Q_2^T \cdots Q_{n-2}^T
$$

If we define 

$$
Q = Q_{n-2} \cdots Q_2 Q_1,
$$ 

then $H = Q A Q^T$, which is the upper Hessenberg matrix we seek.

Note that the standard form is often written as

$$
H = Q^T A Q, \quad \text{or equivalently} \quad
A = Q H Q^T,
$$

where $Q = Q_1^T Q_2 \cdots Q_{n-2}^T$. This is just a matter of convention.

### Computational Cost

* There are $n-2$ main steps in this reduction.
* At each step $k$, applying the Householder reflector $Q_k$ from the left and the right to the $n \times n$ matrix costs $O(n(n-k))$ flops. The dominant cost at each step is $O(n^2)$.
* The total computational cost is the sum of these steps: $\sum_{k=1}^{n-2} O(n(n-k)) \approx O(n^3)$.
(A more careful analysis shows the cost to be $\frac{10}{3}n^3$ flops if $Q$ is explicitly formed, or $\frac{4}{3}n^3$ flops if only $H$ is needed).

This $O(n^3)$ cost is a one-time, upfront investment. Once we have $H$, each subsequent QR iteration (which consists of an $O(n^2)$ QR factorization and an $O(n^2)$ matrix multiplication) is "fast," allowing the overall algorithm to converge efficiently.