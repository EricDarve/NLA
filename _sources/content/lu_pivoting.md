# LU Factorization with Row Pivoting

## The Analysis of Numerical Error: Forward and Backward Perspectives

When we analyze a numerical algorithm, our central concern is the nature of the error. Because we operate in finite-precision arithmetic, the computed solution will almost never be the exact solution. The critical question is how we should quantify this discrepancy. There are two fundamental ways to approach this inquiry: **forward error analysis** and **backward error analysis**.

Let us represent a general numerical problem as the computation of an output $x$ from some input data $d$. We can model this as the application of a function $f$ to the data, such that $x = f(d)$. For the problem of solving a linear system, our data is the pair $(A, b)$, and our solution is $x$, so we may write $x=f(A,b)$.

Due to the accumulation of rounding errors at each step of a computation, any practical algorithm computes not $f$, but an approximation we might call $\tilde{f}$. The result we obtain is the computed solution $\tilde{x} = \tilde{f}(A,b)$.

### Forward Error Analysis: A Direct but Impractical Question

The most natural question one can ask is: **"How close is the computed answer to the true answer?"**

This is the very essence of **forward error analysis**. We seek to establish a bound on the difference between the computed solution $\tilde{x}$ and the true solution $x$. We can measure this **forward error** in either absolute or relative terms:

$$\text{Absolute Forward Error} = \|\tilde{x} - x\| \quad \text{or} \quad \text{Relative Forward Error} = \frac{\|\tilde{x} - x\|}{\|x\|}$$

While this is the quantity we ultimately wish to understand, determining it directly is often impractical. The principal difficulty is that **it requires knowledge of the true solution $x$**. To compute the forward error $\|\tilde{x} - x\|$, one must first find $x$ by solving the original problem $Ax=b$. This is a computationally intensive task, typically requiring $\mathcal{O}(n^3)$ operations. To analyze the error in our solution, we would have to solve the problem again, perhaps with much higher precision, which defeats the purpose of an efficient analysis.

Moreover, a direct forward error analysis would demand that we meticulously track the propagation of every rounding error through the entirety of the algorithm—a hopelessly complex endeavor.

### Backward Error Analysis: The Easy Path of 'Verifying'

Backward error analysis takes a profoundly different and more insightful approach. Instead of asking how incorrect our answer is for the original problem, it asks: **"Is our computed answer the *exact* solution to a slightly different problem?"**

This approach shifts our perspective from the computational difficulty of *solving* to the relative ease of *verifying*. Given a computed solution $\tilde{x}$, it is computationally inexpensive to assess how well it satisfies the original equation by computing the **residual vector** $r$:

$$r = b - A\tilde{x}$$

This requires only a matrix-vector multiplication and a vector subtraction, an $\mathcal{O}(n^2)$ operation—far cheaper than the original $\mathcal{O}(n^3)$ solve. The residual $r$ is of paramount importance because it tells us precisely which nearby problem our algorithm has solved. A simple rearrangement of the equation reveals that:

$$A\tilde{x} = b - r$$

This means our computed solution $\tilde{x}$ is the **exact** solution to the perturbed problem where the right-hand side is $\tilde{b} = b-r$.

The objective of **backward error analysis** is to demonstrate that this perturbation—in this case, the residual $r$—is small relative to the original data. If we can prove that $\|r\|$ is small, we have shown that our algorithm is **backward stable**. A backward stable algorithm produces an answer that is, in a sense, as good as the input data deserves. It has found the exact solution to a problem that differs from the original by a very small amount, and we are able to establish this with a simple verification step, entirely avoiding the need to compute the true solution $x$. This makes backward error analysis the preferred tool for analyzing the stability of algorithms in numerical linear algebra.


## Sensitivity Analysis: Bridging Backward and Forward Error

We have now established a powerful framework: backward error analysis allows us to prove that our algorithm found the *exact* solution to a *nearby* problem. This is an elegant and practical method for assessing an algorithm's stability.

However, our ultimate goal remains to bound the **forward error**. The missing component is an understanding of how errors in the input data affect the output solution. This is the domain of **sensitivity analysis**.

The central question of sensitivity analysis is:

> If we introduce a small perturbation to the input of a problem, how large is the corresponding perturbation in the output?

This is an intrinsic property of the **problem itself**, wholly independent of the algorithm used for its solution. A problem is **sensitive** or **ill-conditioned** if small perturbations in the input can induce large changes in the output. Conversely, a problem is **insensitive** or **well-conditioned** if small input perturbations lead to proportionally small output changes.

### Formalizing Sensitivity

Let us return to our abstract problem $x = f(d)$. Backward error analysis informs us that our computed solution, $\tilde{x}$, is the exact solution for a perturbed input, $\tilde{d}$.

$$\tilde{x} = f(\tilde{d})$$

We have the **backward error**, $\|d - \tilde{d}\|$, which a good algorithm ensures is small. We want to bound the **forward error**, $\|x - \tilde{x}\| = \|f(d) - f(\tilde{d})\|$.

Sensitivity is the factor that connects them. The **absolute sensitivity** can be defined as the limit of the ratio of the change in output to a small change in input:

$$\text{Sensitivity} = \lim_{\varepsilon \to 0} \;\, \sup_{\|\delta d\| \le \varepsilon} \,\frac{\|f(d + \delta d) - f(d)\|}{\|\delta d\|}$$

For small perturbations, this definition gives rise to the fundamental relationship:

$$\|f(d) - f(\tilde{d})\| \approx (\text{Sensitivity}) \times \|d - \tilde{d}\|$$

Or, stated more simply:

$$\text{Forward Error} \approx (\text{Sensitivity}) \times (\text{Backward Error})$$

This equation is the key. It reveals that the forward error is governed by two independent factors:

1.  The **Backward Error**: This is a property of the **algorithm's stability**. A backward stable algorithm yields a small backward error.
2.  The **Sensitivity**: This is a property of the **problem itself**. It acts as an amplification factor.

If a problem has high sensitivity (i.e., it is ill-conditioned), then even a perfectly backward-stable algorithm can produce a solution with a large forward error. The algorithm has performed its task flawlessly—it found an exact solution to a nearby problem—but the nature of the problem itself dictates that even "nearby" is not close enough to guarantee an accurate solution.

## Condition Number: Sensitivity for Linear Systems

With the general framework of sensitivity established, let us apply it to our core problem: solving the linear system $Ax=b$. This will lead us directly to one of the most important concepts in numerical linear algebra: the **condition number** of a matrix.

Our objective is to understand how small perturbations in the input data, $A$ and $b$, affect the solution, $x$.

### Perturbation of the Right-Hand Side `b`

Let us begin with the simplest case. Assume the matrix $A$ is known exactly, but we introduce a small error $\delta b$ into the right-hand side vector $b$. Our original system is $Ax=b$, and the perturbed system is:

$$A(x+\delta x) = b + \delta b$$

where $\delta x$ is the resulting perturbation in the solution. Since $Ax = b$, we may subtract the original equation to find the relationship governing the error:

$$A\delta x = \delta b \quad \implies \quad \delta x = A^{-1}\delta b$$

To measure the magnitude of this error, we take vector norms:

$$\|\delta x\| \le \|A^{-1}\| \|\delta b\|$$

While this absolute error bound is useful, a **relative error** bound is often more insightful. To obtain the relative error in the solution, $\|\delta x\|/\|x\|$, we use the inequality derived from the original system, $\|b\| = \|Ax\| \le \|A\|\|x\|$, which implies $1/\|x\| \le \|A\|/\|b\|$.

Combining these inequalities provides our measure of sensitivity:

$$\frac{\|\delta x\|}{\|x\|} \le \|A^{-1}\| \frac{\|\delta b\|}{\|x\|} \le \|A^{-1}\| \left( \frac{\|A\|}{\|b\|} \right) \|\delta b\| = (\|A\|\|A^{-1}\|) \frac{\|\delta b\|}{\|b\|}$$

This derivation has naturally revealed the amplification factor that relates the relative error in the input $b$ to the relative error in the output $x$. This factor is the **condition number of A**.

### Perturbation of the Matrix `A`

Next, let us consider the case where $b$ is known exactly, but the matrix $A$ contains a small error $\delta A$. The perturbed system is:

$$(A+\delta A)(x+\delta x) = b$$

The analysis is slightly more complex. Assuming the perturbations $\delta A$ and $\delta x$ are small, we can neglect their product $(\delta A)(\delta x)$ as a second-order term. Expanding the equation yields:

$$Ax + A\delta x + (\delta A)x + \mathcal{O}(\|\delta A\|\|\delta x\|) \approx b$$

Since $Ax=b$, this simplifies to $A\delta x \approx -(\delta A)x$, which gives $\delta x \approx -A^{-1}(\delta A)x$. Taking norms, we find:

$$\|\delta x\| \le \|A^{-1}\| \|\delta A\| \|x\|$$

To obtain the relative error, we divide by $\|x\|$:

$$\frac{\|\delta x\|}{\|x\|} \le \|A^{-1}\|\|\delta A\| = (\|A\|\|A^{-1}\|) \frac{\|\delta A\|}{\|A\|}$$

It is a remarkable result that the exact same amplification factor, $\|A\|\|A^{-1}\|$, appears again.

### The Condition Number

These derivations demonstrate that for the problem of solving $Ax=b$, the sensitivity to perturbations in both $A$ and $b$ is governed by the same quantity. We define this as the **condition number of the matrix A**, denoted $\kappa(A)$:

$$\kappa(A) = \|A\| \|A^{-1}\|$$

The condition number is the formal measure of sensitivity for a linear system. It quantifies the maximum extent to which relative errors in the input data can be magnified in the solution.

* If $\kappa(A)$ is small (close to 1), the matrix is **well-conditioned**.
* If $\kappa(A)$ is large, the matrix is **ill-conditioned**.

A fundamental property is that $\kappa(A) \ge 1$ for any matrix $A$ and any induced matrix norm.

Geometrically, the condition number indicates how severely the linear transformation represented by $A$ distorts the space. A well-conditioned matrix transforms a sphere into a mildly eccentric ellipsoid. An ill-conditioned matrix transforms a sphere into a highly elongated, "cigar-shaped" ellipsoid. For such a matrix, a small change in the vector $b$ can result in a large displacement of the solution $x$, as one attempts to pinpoint a location on this extremely stretched object.

This is why the condition number is of such deep importance. It is an intrinsic property of the problem matrix $A$. Regardless of the stability of our algorithm, if $\kappa(A)$ is large, we must anticipate the possibility of an inaccurate solution.


## The General Perturbation Theorem

We have now explored the individual concepts of forward error, backward error, and the condition number. The following theorem synthesizes these ideas into a single, powerful bound that is a cornerstone of numerical linear algebra. It provides a worst-case bound on the forward error of a computed solution as a function of the backward error and the problem's conditioning.

````{prf:theorem} Perturbation Bound for Linear Systems
:label: thm:perturbation_bound
Let $x$ be the solution to $Ax=b$ and $\tilde{x}$ be the solution to the perturbed system $(A+E)\tilde{x} = b+e$. Assume $A$ is invertible and that $\|A^{-1}\|\|E\| < 1$. Then the relative forward error is bounded by:

$$
\frac{\| \tilde{x} - x \|}{\|x\|} \le
\frac{\kappa(A)}{1 - \kappa(A) \frac{\|E\|}{\|A\|}}
\left( \frac{\|E\|}{\|A\|} + \frac{\|e\|}{\|b\|} \right)
$$

where $\kappa(A) = \|A\|\|A^{-1}\|$ is the condition number of the matrix $A$.
````

````{prf:proof}
The proof systematically derives an expression for the error $\delta x = \tilde{x}-x$ and then bounds its norm.

1.  Start with the perturbed system:
    $(A+E)\tilde{x} = b+e$

2.  Isolate terms involving the true matrix A:
    $A\tilde{x} = b+e - E\tilde{x}$

3.  Introduce the true solution by substituting $b = Ax$:
    $A\tilde{x} = Ax + e - E\tilde{x}$

4.  Rearrange to find an expression for the error:
    $A(\tilde{x}-x) = e - E\tilde{x}$.

    Let $\delta x = \tilde{x}-x$. We can write $\tilde{x} = x + \delta x$.
    
    $$A\delta x = e - E(x + \delta x) = e - Ex - E\delta x$$

5.  Solve for the error term $\delta x$: $(A+E)\delta x = e - Ex.$ Hence,

    $$\delta x = (A+E)^{-1} (e - Ex)$$

    This requires $(A+E)$ to be invertible. The condition $\|A^{-1}\|\|E\| < 1$ guarantees this, via the Banach lemma (see note below).

6.  Take norms to bound the error:
    
    $$\|\delta x\| \le \|(A+E)^{-1}\| (\|e\| + \|E\|\|x\|)$$

    We can write 
    
    $$(A+E)^{-1} = (A(I+A^{-1}E))^{-1} = (I+A^{-1}E)^{-1}A^{-1}.$$

    Using the inequality $\|(I+X)^{-1}\| \le \frac{1}{1-\|X\|}$ for $\|X\|<1$ (see note below), we get:
    
    $$\|\delta x\| \le \frac{\|A^{-1}\|}{1-\|A^{-1}E\|} (\|e\| + \|E\|\|x\|)$$

7.  Simplify and introduce relative errors:
    Since $\|A^{-1}E\| \le \|A^{-1}\|\|E\|$, we can weaken the bound slightly to simplify:
    
    $$\|\delta x\| \le \frac{\|A^{-1}\|}{1 - \|A^{-1}\|\|E\|} (\|e\| + \|E\|\|x\|)$$

    Now, divide by $\|x\|$ to get the relative forward error:
    
    $$\frac{\|\delta x\|}{\|x\|} \le \frac{\|A^{-1}\|}{1 - \|A^{-1}\|\|E\|} \left( \frac{\|e\|}{\|x\|} + \|E\| \right)$$

8.  Introduce $\kappa(A)$ and relative backward errors.
    We use the fact that 
    
    $$\|b\| \le \|A\|\|x\| \implies \frac{1}{\|x\|} \le \frac{\|A\|}{\|b\|}.$$

    Thus, we have

    $$\frac{\|\delta x\|}{\|x\|} \le \frac{\|A^{-1}\|}{1 - \|A^{-1}\|\|E\|} \left( \frac{\|e\|}{\|b\|}\|A\| + \|E\| \right)$$

    Finally, multiply the numerator and denominator of the fraction by $\|A\|$:

    $$\frac{\|\delta x\|}{\|x\|} \le \frac{\|A\|\|A^{-1}\|}{1 - \|A\|\|A^{-1}\|\frac{\|E\|}{\|A\|}} \left( \frac{\|e\|}{\|b\|} + \frac{\|E\|}{\|A\|} \right)$$

This completes the proof.
````

### Interpretation

This theorem is the culmination of our analysis. It demonstrates that the relative forward error is bounded by the product of the condition number and the sum of the relative backward errors, adjusted by a factor in the denominator.

If an algorithm is **backward stable**, the relative backward errors $\|E\|/\|A\|$ and $\|e\|/\|b\|$ will be on the order of the machine unit roundoff, $u$. In such cases, the term $\kappa(A)\frac{\|E\|}{\|A\|}$ in the denominator is typically negligible compared to 1, and the bound simplifies to our essential rule of thumb:

$$\text{Relative Forward Error} \lesssim \kappa(A) \times (\text{Relative Backward Error})$$

### Mathematical Aside: The Banach Lemma and Neumann Series

The inequality used in the proof, $\|(I+X)^{-1}\| \le \frac{1}{1-\|X\|}$, is a direct consequence of a result known as the **Banach lemma**. 

````{prf:lemma} Banach Lemma
:label: lem:banach
The Banach lemma states that if $\|X\| < 1$ for some submultiplicative matrix norm, then $(I+X)$ is invertible.
````

````{prf:proof}
The intuition comes from the familiar geometric series for scalars: if $|r|<1$, then 

$$(1+r)^{-1} = 1 - r + r^2 - r^3 + \dots.$$

We can propose a similar series for the matrix inverse, called the Neumann series:

$$(I+X)^{-1} = I - X + X^2 - X^3 + \dots = \sum_{k=0}^{\infty} (-X)^k$$

This series converges precisely because we assume $\|X\| < 1$. Now, we can take the norm of this series to derive the inequality:

1.  Take the norm of the series:
    
    $$\|(I+X)^{-1}\| = \left\| \sum_{k=0}^{\infty} (-X)^k \right\|$$

2.  Apply the triangle inequality: The norm of a sum is less than or equal to the sum of the norms.

    $$\left\| \sum_{k=0}^{\infty} (-X)^k \right\| \le \sum_{k=0}^{\infty} \|(-X)^k\| = \sum_{k=0}^{\infty} \|X\|^k$$

3.  Sum the geometric series: The final term is a standard scalar geometric series with ratio $\|X\| < 1$.

    $$\sum_{k=0}^{\infty} \|X\|^k = 1 + \|X\| + \|X\|^2 + \dots = \frac{1}{1-\|X\|}$$

This gives us the desired result:

$$\|(I+X)^{-1}\| \le \frac{1}{1-\|X\|}$$
````

This lemma is a fundamental tool for perturbation theory, as it allows us to guarantee the invertibility of a perturbed matrix and to bound the norm of its inverse.

## Backward Error Analysis of LU Factorization

Having established a robust framework for analyzing algorithms, we now apply it to our primary algorithm for solving dense linear systems: **LU factorization**. The process of solving $Ax=b$ via LU involves two main stages: first, the factorization $A \approx \tilde{L}\tilde{U}$, and second, the solution of two triangular systems, $\tilde{L}y=b$ and $\tilde{U}\tilde{x}=y$.

A detailed analysis, which we will not reproduce here, consolidates the rounding errors from all of these steps into a single backward error result. The computed solution $\tilde{x}$ is shown to be the exact solution of a nearby system $(A+E)\tilde{x} = b$. The crucial question is: how large is the backward error matrix $E$?

### The Backward Error Bound

````{prf:theorem} Backward Error of LU Factorization
:label: thm:backward_error_lu
The result of a careful backward error analysis of the LU factorization is the following component-wise bound on the error matrix $E$:

$$|E| \le n u (2|A| + 4 |\tilde{L}| |\tilde{U}|) + \mathcal{O}(u^2)$$

Here, $|\cdot|$ denotes the matrix of absolute values of the entries, $u$ is the unit roundoff, and $\tilde{L}$ and $\tilde{U}$ are the computed LU factors of $A$.
````

At first glance, this bound appears complex, but its message is clear. The stability of the algorithm depends entirely on the magnitude of the entries in the computed factors $\tilde{L}$ and $\tilde{U}$. If the entries in $|\tilde{L}|$ and $|\tilde{U}|$ are of the same order of magnitude as the entries in $|A|$, then $|E|$ will be on the order of $u|A|$, and the algorithm is **backward stable**.

However, if significant element growth occurs, and the entries of $|\tilde{L}||\tilde{U}|$ become much larger than those of $|A|$, the backward error $|E|$ can become arbitrarily large. This leads to a critical conclusion: 

:::{admonition} Stability of LU Factorization
:class: important

LU factorization without a proper pivoting strategy is **not** a backward stable algorithm.
:::

### The Mechanism of Instability: Small Pivots

This instability is not merely a theoretical concern. Consider the matrix:

$$
A =
\begin{pmatrix}
\epsilon & 1 \\
1 & \pi
\end{pmatrix}
$$

The exact LU factors are:

$$
L = \begin{pmatrix}
1 & 0 \\ \epsilon^{-1} & 1
\end{pmatrix}, \qquad
U = \begin{pmatrix}
\epsilon & 1 \\ 0 & \pi - \epsilon^{-1}
\end{pmatrix}
$$

For any small $\epsilon > 0$, the entry $\epsilon^{-1}$ in the $L$ factor becomes enormous. When this is substituted into our error bound, the term $|\tilde{L}||\tilde{U}|$ will be dominated by this large element, implying the backward error $\|E\|$ will be on the order of $u \cdot \epsilon^{-1}$. By choosing $\epsilon$ to be very small (e.g., on the order of machine precision), we can make the backward error arbitrarily large. An algorithm that produces an $\mathcal{O}(1)$ backward error for a well-conditioned problem is numerically useless.

The source of these large entries is evident from the mechanics of Gaussian elimination. At each step $k$, we compute multipliers by dividing by the pivot element $a_{kk}$:

$$l_{ik} = \frac{a_{ik}}{a_{kk}}, \quad \text{for } i > k$$

If a pivot $a_{kk}$ is very small compared to the entries below it in the same column, the resulting multipliers in $L$ will be very large, leading directly to the instability captured in the error bound. This is precisely the issue that pivoting strategies are designed to prevent.

### An Intuitive Parallel with Summation Error

This result for the backward error of LU factorization should feel familiar. It's a direct reflection of our earlier analysis of roundoff errors in simple summation, where we found the error bound (see {prf:ref}`thm:summation_error`):

$$|\Delta S_n| \le n u \sum_{i=1}^n |x_i| + \mathcal{O}(u^2)$$

This bound taught us a crucial lesson: the accuracy of a sum depends on the **magnitude of the numbers being added**, not just the magnitude of the final result.

We can see the exact same principle at play in the LU factorization. The process of factorization and, conversely, the reconstruction of the matrix $A$ from its factors, relies on inner products:

$$a_{ij} = \sum_{k=1}^n l_{ik} u_{kj}$$

Each entry of $A$ is the result of a sum. If the entries of $L$ and $U$ are large, then the individual terms $l_{ik}u_{kj}$ in this sum will be large. Our summation error bound tells us that computing this sum will likely incur a large absolute error.

The term $|\tilde{L}||\tilde{U}|$ in the backward error bound for $E$ is precisely the matrix-level consequence of this fundamental observation. Its entries, 

$$(|\tilde{L}||\tilde{U}|)_{ij} = \sum_k |\tilde{l}_{ik}||\tilde{u}_{kj}|,$$

represent the worst-case sum of magnitudes for the terms that form $a_{ij}$. Therefore, the instability we observe when element growth occurs is perfectly consistent with our understanding of basic floating-point operations. Large intermediate numbers are a universal sign of potential numerical instability.

## The Solution: LU Factorization with Partial Pivoting

Now that we have identified the cause of instability in the basic LU algorithm, we can introduce the standard remedy.

### Partial Pivoting

The instability in Gaussian elimination arises when a pivot element $a_{kk}$ is small relative to the entries below it. The solution is both simple and elegant: at each step, we must ensure that the pivot is as large as possible.

This strategy is called **partial pivoting** (or row pivoting). At each step $k$ of the elimination, the algorithm is modified as follows:

1.  **Search for the largest pivot:** Before performing elimination for column $k$, find the entry in that column with the largest absolute value, on or below the diagonal. Let this be $a_{pk}$, where $p \ge k$.
2.  **Swap rows:** Interchange row $p$ with the current pivot row, row $k$.
3.  **Eliminate:** Proceed with elimination as usual, now using the largest possible pivot for that column.

By construction, this strategy guarantees that the pivot $a_{kk}$ is the largest entry (in magnitude) in its column among all rows not yet used as pivot rows. This has a crucial consequence for the multipliers:

$$|l_{ik}| = \frac{|a_{ik}|}{|a_{kk}|} \le 1 \quad \text{for } i > k$$

This simple change prevents the uncontrolled growth of elements in the $L$ factor.

### The $PA = LU$ Factorization

The systematic swapping of rows can be represented mathematically by a **permutation matrix** $P$. A permutation matrix is an identity matrix with its rows reordered. Pre-multiplying a matrix $A$ by $P$ (i.e., forming $PA$) has the effect of applying those same row permutations to $A$.

Therefore, LU factorization with partial pivoting does not compute the factors of $A$ itself, but rather of a row-permuted version of $A$. The resulting factorization is:

$$PA = LU$$

This is the standard form of LU factorization implemented in virtually all numerical software. It comes with a powerful guarantee that the basic version lacks: the $PA=LU$ factorization **exists for any square matrix A**, singular or not (see below for a proof).

### Implementing LU with Partial Pivoting

Below is a simple implementation of LU factorization with partial pivoting in Python. The function modifies the input matrix `A` in place to store the factors `L` and `U`, and returns the permutation matrix `P`.

```python
import numpy as np

def lu_factorization_with_row_pivoting(A):
    """
    Perform LU factorization with row (partial) pivoting in place.
    On exit, A stores L (unit lower, diagonal = 1 implicit) in its strictly 
    lower triangle and U in its upper triangle. 
    Returns P such that P @ A_orig = L @ U.
    """
    n = A.shape[0]
    P = np.eye(n, dtype=A.dtype)

    for k in range(n - 1):
        # Pivot: index of max |A[i,k]| for i >= k
        p = k + np.argmax(np.abs(A[k:, k]))
        if p != k:
            A[[k, p], :] = A[[p, k], :]
            P[[k, p], :] = P[[p, k], :]

        if A[k, k] == 0:
            continue # Skip elimination

        # Update the k-th column of L
        A[k+1:n, k] /= A[k, k]
        # Rank-one update of the trailing submatrix
        A[k+1:n, k+1:n] -= np.outer(A[k+1:n, k], A[k, k+1:n])

    return P
```

### The Unstable Example, Stabilized

Let us apply this strategy to the matrix that previously demonstrated instability:

$$
A =
\begin{pmatrix}
\epsilon & 1 \\
1 & \pi
\end{pmatrix}
$$

At the first step, the pivot candidates in the first column are $\epsilon$ and $1$. Since $|1| > |\epsilon|$, the partial pivoting strategy mandates a swap of row 1 and row 2. The permutation matrix and the resulting permuted matrix are:

$$P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad PA = \begin{pmatrix} 1 & \pi \\ \epsilon & 1 \end{pmatrix}$$

Now, we perform LU factorization on $PA$:

$$
L = \begin{pmatrix}
1 & 0 \\ \epsilon & 1
\end{pmatrix}, \qquad
U = \begin{pmatrix}
1 & \pi \\ 0 & 1 - \epsilon \pi
\end{pmatrix}
$$

The enormous $\epsilon^{-1}$ term has vanished. All entries in $L$ and $U$ are of moderate size. The backward error is now proportional to $u$, rendering the algorithm **backward stable** for this problem.

### Proof of the Existence of LU Factorization with Row Pivoting

The use of partial pivoting not only stabilizes the LU factorization algorithm against the growth of roundoff error, but it also provides a powerful theoretical guarantee: the factorization is guaranteed to exist for *any* square matrix, whether it is invertible or not. This is a significant improvement over the basic LU factorization, which can fail even for invertible matrices if a zero pivot is encountered.

```{prf:theorem} Existence of LU Factorization with Row Pivoting
:label: thm:existence_lu_pivoting
For any square matrix $A\in\mathbb{F}^{n\times n}$ (where $\mathbb{F}=\mathbb{R}$ or $\mathbb{C}$), there exist a permutation matrix $P$, a unit lower–triangular matrix $L$, and an upper–triangular matrix $U$ such that

$$
PA = LU .
$$
```

```{prf:proof}
We proceed by induction on the dimension $n$. For notational simplicity, we use $\mathbb F = \mathbb R$; the complex case is identical.

**Base case ($n=1$):**
For a $1 \times 1$ matrix $A=[a]$, the factorization is trivial. We can choose $P=[1]$, $L=[1]$, and $U=[a]$. Then $PA = [1][a] = [a]$ and $LU = [1][a] = [a]$, so the statement holds.

**Inductive step:**
Assume the statement holds for all matrices of size $(n-1)\times(n-1)$. We must show it holds for an arbitrary $n\times n$ matrix $A$.

1.  **Choose a pivot and permute.**
    Search the first column of $A$ for an element with the largest absolute value. Let this element be in row $p$. If the entire first column is zero, we can choose any row (e.g., $p=1$). Let $P_1$ be the permutation matrix that swaps row 1 and row $p$. We then form the permuted matrix:

    $$
    A^{(1)} := P_1 A =
    \begin{bmatrix}
    \alpha & w^T \\
    x & A_{22}
    \end{bmatrix}
    $$

    where $\alpha \in \mathbb{R}$ is the pivot element, $w, x \in \mathbb{R}^{n-1}$, and $A_{22} \in \mathbb{R}^{(n-1)\times(n-1)}$. By our choice of pivot, $|\alpha| \ge |x_i|$ for all entries $x_i$ in the vector $x$.

2.  **Eliminate below the pivot.**
    We now perform one step of elimination. Define the vector of multipliers $\ell$ as:

    $$
    \ell :=
    \begin{cases}
    \alpha^{-1}x, & \text{if }\alpha\neq 0,\\
    0,            & \text{if }\alpha=0.
    \end{cases}
    $$

    Note that if $\alpha=0$, our pivot choice implies that the entire first column is zero, so $x=0$, and $\ell=0$ is the correct choice. We can write this elimination step using a matrix $E_1$:

    $$
    E_1 =
    \begin{bmatrix}
    1 & 0 \\
    -\ell & I
    \end{bmatrix}
    \quad \implies \quad
    E_1 A^{(1)} =
    \begin{bmatrix}
    \alpha & w^T \\
    0 & S
    \end{bmatrix}
    $$

    where $S = A_{22} - \ell w^T$ is the Schur complement. The matrix $E_1$ is unit lower-triangular, and its inverse $E_1^{-1} = \begin{bmatrix}1&0\\ \ell & I\end{bmatrix}$ is also unit lower-triangular.

3.  **Apply the induction hypothesis.**
    The Schur complement $S$ is an $(n-1)\times(n-1)$ matrix. By our induction hypothesis, there exists a factorization for $S$:

    $$
    P_S S = L_S U_S
    $$

    where $P_S$ is an $(n-1)\times(n-1)$ permutation matrix, $L_S$ is unit lower-triangular, and $U_S$ is upper-triangular.

4.  **Assemble the final factorization.**
    We can "lift" the permutation $P_S$ to the full $n \times n$ dimension by defining a block matrix $\Pi$:

    $$
    \Pi := \begin{bmatrix} 1 & 0 \\ 0 & P_S \end{bmatrix}
    $$

    Now, we multiply our previous result by $\Pi$:

    $$
    \Pi (E_1 A^{(1)}) =
    \begin{bmatrix} 1 & 0 \\ 0 & P_S \end{bmatrix}
    \begin{bmatrix} \alpha & w^T \\ 0 & S \end{bmatrix}
    =
    \begin{bmatrix} \alpha & w^T \\ 0 & P_S S \end{bmatrix}
    $$

    Substituting the factorization for $P_S S$, we get:

    $$
    \Pi E_1 P_1 A =
    \begin{bmatrix} \alpha & w^T \\ 0 & L_S U_S \end{bmatrix}
    =
    \underbrace{\begin{bmatrix} 1 & 0 \\ 0 & L_S \end{bmatrix}}_{\text{unit lower-triangular}} \hspace{1em}
    \underbrace{\begin{bmatrix} \alpha & w^T \\ 0 & U_S \end{bmatrix}}_{\text{upper-triangular}}
    $$

    Let's call these two new matrices $\tilde{L}$ and $U$. We now have $\Pi E_1 P_1 A = \tilde{L} U$. To isolate $A$, we rearrange:

    $$
    (\Pi P_1) A = (\Pi E_1^{-1} \Pi^{-1}) \tilde{L} U
    $$

    Let's examine the terms:
    -   $P = \Pi P_1$ is a product of permutation matrices, and is therefore a permutation matrix.
    -   $L = (\Pi E_1^{-1} \Pi^{-1}) \tilde{L}$. The matrix $E_1^{-1} = \begin{bmatrix}1&0\\ \ell & I\end{bmatrix}$ is unit lower-triangular. Conjugating by $\Pi$ only permutes rows and columns 2 through $n$, which preserves the unit lower-triangular structure. The product of two unit lower-triangular matrices ($\Pi E_1^{-1} \Pi^{-1}$ and $\tilde{L}$) is also unit lower-triangular.
    -   $U$ is already upper-triangular by construction.

    We have successfully constructed $P, L, U$ with the required properties such that $PA=LU$. This completes the induction.
```

**Remarks:**
- This proof is constructive and directly corresponds to the LU factorization algorithm with row pivoting.
- The proof does not require the matrix $A$ to be invertible. If $A$ is singular, the process still works, and the singularity will manifest as one or more zero entries on the diagonal of the upper-triangular factor $U$.

### A Final Word on Stability

It is crucial to note that while partial pivoting guarantees $|l_{ik}| \le 1$, it does *not* offer a mathematical guarantee that the entries of the $U$ factor will not grow large. It is possible to construct matrices where significant element growth still occurs in $U$.

An example of a matrix that can cause large element growth is provided below, along with a discussion of the theoretical worst-case growth factor.

However, decades of practical experience have shown that such cases are exceptionally rare. For the vast majority of problems encountered in science and engineering, LU with partial pivoting is a remarkably robust and accurate algorithm. It is the gold standard for solving dense linear systems directly.

### An Example of Element Growth

A classic example of a matrix that exhibits significant element growth in its $U$ factor, even with partial pivoting, is one with 1s on the diagonal, -1s in the lower triangle, and 1s in the last column.

The theoretical worst-case growth factor for an $n \times n$ matrix is $2^{n-1}$.

Consider the following $4 \times 4$ matrix:

$$
A =
\begin{pmatrix}
\phantom{-}1 & \phantom{-}0 & \phantom{-}0 & 1 \\
-1 & \phantom{-}1 & \phantom{-}0 & 1 \\
-1 & -1 & \phantom{-}1 & 1 \\
-1 & -1 & -1 & 1
\end{pmatrix}
$$

When you apply LU factorization with **row pivoting** to this matrix, the pivoting strategy does nothing. At each step, the largest element in the column is already on the diagonal (its magnitude is 1), so no row swaps occur.

However, the elimination steps cause the entries in the last column to double at each stage. The resulting upper triangular matrix $U$ is:

$$
U =
\begin{pmatrix}
1 & 0 & 0 & 1 \\
0 & 1 & 0 & 2 \\
0 & 0 & 1 & 4 \\
0 & 0 & 0 & 8
\end{pmatrix}
$$

with

$$ L = \begin{pmatrix} \phantom{-}1 & \phantom{-}0 & \phantom{-}0 & 0 \\ -1 & \phantom{-}1 & \phantom{-}0 & 0 \\ -1 & -1 & \phantom{-}1 & 0 \\ -1 & -1 & -1 & 1 \end{pmatrix}$$

The largest element in the original matrix $A$ was 1, but the largest element in the computed factor $U$ has grown to 8. This demonstrates significant element growth despite the use of row pivoting.

### The Worst-Case Growth Factor

The example above illustrates the general worst-case scenario. For an $n \times n$ matrix of this type, the bottom-right element of the matrix $U$ will grow to $2^{n-1}$.

The **growth factor** is defined as the ratio of the largest element (in magnitude) in the computed factors to the largest element in the original matrix:

$$\rho_n = \frac{\max_{i,j} |u_{ij}|}{\max_{i,j} |a_{ij}|}$$

For LU with row pivoting, the proven worst-case bound for this growth factor is:

$$\rho_n = 2^{n-1}$$

This exponential growth is a sobering theoretical result. It implies that for a large matrix (e.g., $n=50$), the entries could theoretically become larger by a factor of $2^{49}$, leading to a catastrophic loss of accuracy.

Fortunately, this worst-case behavior is almost never seen in practice. Matrices that exhibit this exponential growth are highly contrived and do not typically arise from real-world problems. For most practical applications, the growth factor remains small, making LU with row pivoting a very reliable and stable algorithm.