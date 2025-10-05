# Pivoting in LU Factorization

## Forward and Backward Error Analysis: A Tale of Two Questions

When we design and analyze numerical algorithms, our primary concern is understanding the error. Since we are working with finite-precision floating-point arithmetic, our computed solution will almost never be the exact solution. The question is, how do we quantify this error? There are two fundamental ways to approach this, known as **forward error analysis** and **backward error analysis**.

Let's represent a general problem as computing an output $x$ from some input data $d$. We can think of this as applying a function $f$ to the data: $x = f(d)$. For the problem of solving a linear system, our data is the pair $(A, b)$ and our solution is $x$, so $x=f(A,b)$.

Due to rounding errors at each step of our computation, our algorithm doesn't compute $f$, but rather an approximation, which we'll call $\tilde{f}$. The result we actually get is the computed solution $\tilde{x} = \tilde{f}(A,b)$.

### Forward Error Analysis: The Difficult Path of 'Solving'

The most natural question to ask is: **"How close is our computed answer to the true answer?"**

This is the essence of **forward error analysis**. We want to find a bound on the difference between the computed solution $\tilde{x}$ and the true solution $x$. Mathematically, we seek a bound on the **forward error**, which can be measured in absolute or relative terms:

$$\text{Absolute Forward Error} = \|\tilde{x} - x\| \quad \text{or} \quad \text{Relative Forward Error} = \frac{\|\tilde{x} - x\|}{\|x\|}$$

While this is the question we ultimately care about, answering it directly is often impractical because **it requires you to know the true solution $x$**. To compute the forward error $\|\tilde{x} - x\|$, you must first find $x$ by solving the original problem $Ax=b$. This is a computationally intensive task, typically requiring $\mathcal{O}(n^3)$ operations. In essence, to analyze the error of your solution, you have to go through the hard work of solving the entire problem again, for example with higher precision, which defeats the purpose of an efficient analysis.

Furthermore, a direct forward error analysis would require us to meticulously track how every single rounding error propagates through the entire algorithm. This is a hopelessly complex task.

### Backward Error Analysis: The Easy Path of 'Verifying'

Backward error analysis takes a brilliantly different approach. Instead of asking how wrong our answer is, it asks: **"Is our computed answer the *exact* answer to a slightly different problem?"**

This approach shifts our perspective from the difficulty of *solving* to the ease of *verifying*. Given a computed solution $\tilde{x}$, it is computationally cheap to check how well it satisfies the original equation. We can do this by computing the **residual vector** $r$:

$$r = b - A\tilde{x}$$

This is just a matrix-vector multiplication and a vector subtraction, an $\mathcal{O}(n^2)$ operation—much cheaper than the original $\mathcal{O}(n^3)$ solve. The residual $r$ is crucial because it directly tells us which nearby problem our algorithm *actually* solved. A simple rearrangement of the equation shows:

$$A\tilde{x} = b - r$$

This means our computed solution $\tilde{x}$ is the **exact** solution to the slightly perturbed problem where the right-hand side is $\tilde{b} = b-r$.

The goal of **backward error analysis** is then to show that this perturbation (in this case, the residual $r$) is small relative to the original data. If we can prove that $\|r\|$ is small, we have shown that our algorithm is **backward stable**. A backward stable algorithm delivers an answer that is, in a sense, as good as the input data deserves. It has found the exact solution to a problem that differs from the original by a very small amount, and we were able to prove this with a simple verification step, avoiding the need to solve for the true $x$ at all. This makes backward error analysis the preferred tool for analyzing the stability of algorithms in numerical linear algebra.

## Sensitivity Analysis: The Bridge from Backward to Forward Error

We've established a powerful framework: backward error analysis allows us to prove that our algorithm found the *exact* solution to a *nearby* problem. This is an elegant and practical way to assess an algorithm's stability.

However, our ultimate goal remains to bound the **forward error**—the error in the final answer. The missing piece of the puzzle is understanding how errors in the input data affect the output solution. This is the domain of **sensitivity analysis**.

The central question of sensitivity analysis is:

> If we make a small change to the input of a problem, how much does the output change?

This property is intrinsic to the **problem itself**, regardless of the algorithm used to solve it. A problem is **sensitive** or **ill-conditioned** if small perturbations in the input can cause large changes in the output. A problem is **insensitive** or **well-conditioned** if small input perturbations lead to proportionally small output changes.

### Formalizing Sensitivity

Let's return to our abstract problem $x = f(d)$, where $d$ is the input data and $x$ is the solution. Backward error analysis tells us that our computed solution, $\tilde{x}$, is the exact solution for a perturbed input, $\tilde{d}$.

$$\tilde{x} = f(\tilde{d})$$

We have the **backward error**, which we can bound: $\|d - \tilde{d}\|$.
We want the **forward error**, which we need to find: $\|x - \tilde{x}\| = \|f(d) - f(\tilde{d})\|$.

Sensitivity is the factor that connects them. The **absolute sensitivity** is defined as the maximum ratio of the change in output to a small change in input:

$$\text{Sensitivity} = \lim_{\varepsilon \to 0} \;\, \sup_{\|\delta d\| \le \varepsilon} \,\frac{\|f(d + \delta d) - f(d)\|}{\|\delta d\|}$$

For small perturbations, this gives us the fundamental relationship:

$$\|f(d) - f(\tilde{d})\| \approx (\text{Sensitivity}) \times \|d - \tilde{d}\|$$

Or, more simply:

$$\text{Forward Error} \approx (\text{Sensitivity}) \times (\text{Backward Error})$$

This equation is the key. It shows that the forward error we care about is controlled by two independent factors:

1.  The **Backward Error**: This depends on the **stability of our algorithm**. A good algorithm gives a small backward error.
2.  The **Sensitivity**: This depends only on the **problem itself**. It acts as an amplification factor.

If a problem has high sensitivity (it's ill-conditioned), then even a perfectly backward-stable algorithm can produce a solution with a large forward error. The algorithm has done its job perfectly—it found an exact solution to a problem very close by—but the nature of the problem itself means that even "very close" is not close enough to guarantee an accurate answer.

## Condition Number: Sensitivity for Linear Systems

### Perturbing the Right-Hand Side `b`

Excellent. Now that we have the general framework of sensitivity, let's apply it to our core problem: solving the linear system $Ax=b$. This will lead us directly to one of the most important concepts in numerical linear algebra: the **condition number** of a matrix.

Our goal is to understand how small perturbations in the input data, $A$ and $b$, affect the solution, $x$.

Let's start with the simplest case. Assume the matrix $A$ is known exactly, but we introduce a small error $\delta b$ into the right-hand side vector $b$. Our original system is $Ax=b$. The perturbed system is:

$$A(x+\delta x) = b + \delta b$$

where $\delta x$ is the resulting error in the solution. Since $Ax = b$, we can subtract the original equation to find the error relationship:

$$A\delta x = \delta b \quad \implies \quad \delta x = A^{-1}\delta b$$

To measure the size of this error, we take vector norms:

$$\|\delta x\| \le \|A^{-1}\| \|\delta b\|$$

This absolute error bound is useful, but a **relative error** is often more insightful. To get the relative error in the solution, $\|\delta x\|/\|x\|$, we can use the fact that $\|b\| = \|Ax\| \le \|A\|\|x\|$, which implies $1/\|x\| \le \|A\|/\|b\|$.

Combining these inequalities gives us our sensitivity measure:

$$\frac{\|\delta x\|}{\|x\|} \le \|A^{-1}\| \frac{\|\delta b\|}{\|x\|} \le \|A^{-1}\| \left( \frac{\|A\|}{\|b\|} \right) \|\delta b\| = (\|A\|\|A^{-1}\|) \frac{\|\delta b\|}{\|b\|}$$

This derivation naturally reveals the amplification factor that relates the relative error in the input $b$ to the relative error in the output $x$. This factor is the **condition number of A**.

### Perturbing the Matrix `A`

Now, let's consider the case where $b$ is known exactly, but the matrix $A$ has a small error $\delta A$. The perturbed system is:

$$(A+\delta A)(x+\delta x) = b$$

This is slightly more complex to analyze. If we assume the perturbations $\delta A$ and $\delta x$ are small, we can ignore their product $(\delta A)(\delta x)$ as a negligible second-order term. Expanding the equation gives:

$$Ax + A\delta x + (\delta A)x + \mathcal{O}(\|\delta A\|\|\delta x\|) \approx b$$

Since $Ax=b$, this simplifies to $A\delta x \approx -(\delta A)x$, which means $\delta x \approx -A^{-1}(\delta A)x$. Taking norms, we get:

$$\|\delta x\| \le \|A^{-1}\| \|\delta A\| \|x\|$$

To get the relative error, we divide by $\|x\|$:

$$\frac{\|\delta x\|}{\|x\|} \le \|A^{-1}\|\|\delta A\| = (\|A\|\|A^{-1}\|) \frac{\|\delta A\|}{\|A\|}$$

Remarkably, the exact same amplification factor, $\|A\|\|A^{-1}\|$, appears again.

## The Condition Number of a Matrix 👑

These derivations show that for the problem of solving $Ax=b$, the sensitivity to perturbations in both $A$ and $b$ is governed by the same quantity. We define this as the **condition number of the matrix A**, denoted $\kappa(A)$:

$$\kappa(A) = \|A\| \|A^{-1}\|$$

The condition number is the formal measure of sensitivity for solving a linear system. It tells us the maximum extent to which relative errors in the input data can be magnified in the solution.

* If $\kappa(A)$ is small (close to 1), the matrix is **well-conditioned**.
* If $\kappa(A)$ is large, the matrix is **ill-conditioned**.

A key property is that $\kappa(A) \ge 1$ for any matrix.

Geometrically, the condition number tells you how much the matrix $A$ distorts the geometry of the space. A well-conditioned matrix transforms a sphere into a slightly stretched sphere (an ellipsoid with similar axes). An ill-conditioned matrix transforms a sphere into a very long, thin "cigar" shape. For an ill-conditioned matrix, a tiny change in the vector $b$ can result in a huge jump in the solution $x$, because you're trying to pinpoint a location on an extremely elongated ellipse.

This is why we care so deeply about the condition number. It's an intrinsic property of the problem matrix $A$. No matter how stable our algorithm is (i.e., how small we make the backward errors), if $\kappa(A)$ is large, we must be prepared for the possibility of an inaccurate solution.

### The General Perturbation Theorem

We have now explored the individual concepts of forward error, backward error, and the condition number. The following theorem combines them into a single, powerful bound that is central to numerical linear algebra. It provides a worst-case bound on the forward error of a computed solution as a function of the backward error and the problem's condition number.

### Theorem

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

This result is the punchline of our analysis. It shows that the relative forward error is bounded by the product of the condition number and the sum of the relative backward errors, with a small correction factor in the denominator.

If an algorithm is **backward stable**, the relative backward errors $\|E\|/\|A\|$ and $\|e\|/\|b\|$ will be on the order of the unit roundoff, $u$. In this case, the term $\kappa(A)\frac{\|E\|}{\|A\|}$ in the denominator is usually very small compared to 1, and the bound simplifies to our rule of thumb:

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

We have established a robust framework for analyzing algorithms using backward error and conditioning. We now apply this framework to our primary algorithm for solving dense linear systems: **LU factorization**.

The process of solving $Ax=b$ using LU involves two steps:

1.  Factor $A \approx \tilde{L}\tilde{U}$.
2.  Solve the triangular systems $\tilde{L}y=b$ and $\tilde{U}\tilde{x}=y$.

The rounding errors from all of these steps can be consolidated into a single, elegant backward error result. A detailed analysis, which we will not reproduce here, shows that the computed solution $\tilde{x}$ is the exact solution to a nearby system $(A+E)\tilde{x} = b$. The crucial question is: how large is this backward error matrix $E$?

### The Backward Error Bound

````{prf:theorem} Backward Error of LU Factorization
:label: thm:backward_error_lu
The result of a careful backward error analysis of the LU factorization is the following component-wise bound on the error matrix $E$:

$$|E| \le n u (2|A| + 4 |\tilde{L}| |\tilde{U}|) + \mathcal{O}(u^2)$$

Here, $|\cdot|$ denotes the matrix of absolute values of the entries, $u$ is the unit roundoff, and $\tilde{L}$ and $\tilde{U}$ are the computed LU factors of $A$.
````

At first glance, this bound seems complex, but its message is clear. The stability of the algorithm depends entirely on the size of the entries in the computed factors $\tilde{L}$ and $\tilde{U}$. If the entries in $|\tilde{L}|$ and $|\tilde{U}|$ are of the same order of magnitude as the entries in $|A|$, then $|E|$ will be on the order of $u|A|$, and the algorithm is **backward stable**.

However, if element growth occurs and the entries of $|\tilde{L}||\tilde{U}|$ become much larger than $|A|$, the backward error $|E|$ can become arbitrarily large. This leads to a critical conclusion:

```{important} Stability of LU Factorization
LU factorization without a proper pivoting strategy is **not** a backward stable algorithm.
```

### Why Instability Occurs: The Peril of Small Pivots ☢️

This instability is not just a theoretical concern. Let's revisit our previous example:

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

For any small $\epsilon > 0$, the entry $\epsilon^{-1}$ in the $L$ factor becomes enormous. Plugging this into our error bound, the term $|\tilde{L}||\tilde{U}|$ will contain this large element, meaning the backward error $\|E\|$ will be on the order of $u \cdot \epsilon^{-1}$. By choosing $\epsilon$ to be very small (e.g., $\epsilon \approx u$), we can make the backward error very large, perhaps close to $\mathcal{O}(1)$. An algorithm that produces an $\mathcal{O}(1)$ backward error for a well-conditioned problem is, for all practical purposes, useless.

The source of these large entries is clear from the mechanics of the algorithm. At each step $k$, we compute the multipliers in the $k$-th column of $L$ by dividing by the pivot element $a_{kk}$:

$$l_{ik} = \frac{a_{ik}}{a_{kk}}, \quad \text{for } i > k$$

If a pivot $a_{kk}$ happens to be very small compared to the entries below it, the resulting values in $L$ will be very large, leading directly to the instability shown in the error bound. This is precisely the issue that pivoting strategies are designed to prevent.

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

The instability in Gaussian elimination arises when a pivot element $a_{kk}$ is small relative to the entries below it, leading to large multipliers $|l_{ik}| > 1$. The solution is simple and elegant: ensure the pivots are never small.

The strategy is called **partial pivoting** (or row pivoting). At each step $k$ of the elimination process, we modify the algorithm as follows:

1.  **Find the largest pivot:** Before performing the elimination for column $k$, find the entry in that column with the largest absolute value, on or below the diagonal. Let this be $a_{pk}$, where $p \ge k$.
2.  **Swap rows:** Interchange row $p$ with the current pivot row, row $k$.
3.  **Eliminate:** Proceed with the elimination as usual, using the new, largest-possible pivot.

By construction, this strategy guarantees that the pivot $a_{kk}$ is the largest entry (in magnitude) in its column among all rows that have not yet been used as pivot rows. This has a crucial consequence for the multipliers:

$$|l_{ik}| = \frac{|a_{ik}|}{|a_{kk}|} \le 1 \quad \text{for } i > k$$

This simple maneuver completely prevents the uncontrolled growth of elements in the $L$ factor.

### The $PA = LU$ Factorization

The process of systematically swapping rows can be represented mathematically by a **permutation matrix** $P$. A permutation matrix is an identity matrix with its rows reordered. Multiplying a matrix $A$ on the left by $P$ (i.e., $PA$) has the effect of applying the same row permutations to $A$.

Therefore, LU factorization with partial pivoting does not compute the factors of $A$ itself, but rather of a row-permuted version of $A$. The resulting factorization is:

$$PA = LU$$

This is the standard form of LU factorization implemented in virtually all numerical software. It comes with a powerful guarantee that the basic version lacks: the $PA=LU$ factorization **exists for any square matrix A**.

### The Unstable Example, Stabilized ✅

Let's apply this strategy to the matrix that previously gave us trouble:

$$
A =
\begin{pmatrix}
\epsilon & 1 \\
1 & \pi
\end{pmatrix}
$$

At the first step, the pivot candidates are $\epsilon$ and $1$. Since $|1| > |\epsilon|$, the partial pivoting strategy requires us to swap row 1 and row 2. The permutation matrix and the resulting permuted matrix are:

$$
P = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \qquad PA = \begin{pmatrix} 1 & \pi \\ \epsilon & 1 \end{pmatrix}
$$

Now, we perform LU factorization on $PA$:

$$
L = \begin{pmatrix}
1 & 0 \\ \epsilon & 1
\end{pmatrix}, \qquad
U = \begin{pmatrix}
1 & \pi \\ 0 & 1 - \epsilon \pi
\end{pmatrix}
$$

The enormous $\epsilon^{-1}$ term has vanished. All entries in $L$ and $U$ are of moderate size. The backward error is now proportional to $u$, making the algorithm **backward stable** for this problem.

### A Final Word on Stability

It is crucial to note that while partial pivoting guarantees $|l_{ik}| \le 1$, it does *not* offer a mathematical guarantee that the entries of the $U$ factor will not grow large. It is possible to construct matrices where significant element growth still occurs in $U$.

An example of a matrix that can cause large element growth is provided below, along with a discussion of the theoretical worst-case growth factor.

However, decades of practical experience have shown that such cases are exceptionally rare. For the vast majority of problems encountered in science and engineering, LU with partial pivoting is a remarkably robust and accurate algorithm. It is the gold standard for solving dense linear systems directly.

### An Example of Element Growth 📈

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