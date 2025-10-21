# SVD for Rank-Deficient Least-Squares

In our previous discussions, we assumed $A$ had full column rank. But what happens if it doesn't?

* If rank$(A) = r < n$, the matrix $A$ is **rank-deficient**.
* This implies $A$ has a non-trivial **null space**, $N(A)$.
* If $x_p$ is a solution to $\text{argmin}_x \|Ax - b\|_2$, then $x_p + z$ is *also* a solution for *any* vector $z \in N(A)$, because $A(x_p + z) - b = (Ax_p - b) + Az = (Ax_p - b)$.
* This means the least-squares solution is **not unique**.

Our previous methods fail. For example, in the QR factorization $A=QR$, the upper-triangular matrix $R$ will have zeros on its diagonal, making it singular and not invertible.

## The Minimum Norm Solution

When the solution is not unique, we must add a second criterion to select one of the many possible solutions. The standard choice is the **minimum norm solution**, which is the unique vector $x^*$ that satisfies two conditions:

1.  **It solves the LS problem:** $x^* = \text{argmin}_x \|Ax - b\|_2$.
2.  **It has the minimum 2-norm:** $\|x^*\|_2 \le \|x\|_2$ for any other LS solution $x$.

This second condition is equivalent to requiring the solution to be orthogonal to the null space: $x^* \perp N(A)$. This unique solution can be found robustly using the Singular Value Decomposition (SVD).

````{prf:proof}
We prove the fact that the minimum-norm least-squares solution $x^*$ is orthogonal to the null space, $x^* \perp N(A)$.

Let $S$ be the set of all vectors $x$ that solve the least-squares problem. We are looking for a unique vector $x^* \in S$ such that $\|x^*\|_2$ is minimized.

1.  **Decompose any LS solution:**
    Let $x_p$ be *any* solution to the least-squares problem. As we've established, the set of all solutions $S$ can be written as:

    $$
    S = \{ x_p + z \mid z \in N(A) \}
    $$

    where $N(A)$ is the null space of $A$. Our goal is to find the vector in this set with the smallest 2-norm.

2.  **Decompose $x_p$ using Fundamental Subspaces:**
    By the fundamental theorem of linear algebra, any vector in $\mathbb{R}^n$ (like $x_p$) can be uniquely decomposed into a component in the **row space**, $\text{span}(A^T)$, and a component in the **null space**, $N(A)$.

    $$
    x_p = x_r + x_n
    $$

    where $x_r \in \text{span}(A^T)$ and $x_n \in N(A)$.

3.  **Rewrite the set of all solutions:**
    Let's substitute this decomposition back into our set $S$:

    $$
    S = \{ x_r + z \mid z \in N(A) \}
    $$

    since $x_n$ is also in $N(A)$, we can absorb it into the general $z$ term.

4.  **Minimize the Norm using Orthogonality:**
    We want to find $\text{argmin}_{x \in S} \|x\|_2^2$, which is now:

    $$
    \min_{z' \in N(A)} \|x_r + z'\|_2^2
    $$

    A key part of the fundamental theorem is that the row space and the null space are **orthogonal complements**. This means that $x_r \perp z'$ for any $z' \in N(A)$.
    By the Pythagorean theorem, the squared norm becomes:

    $$
    \|x_r + z'\|_2^2 = \|x_r\|_2^2 + \|z'\|_2^2
    $$

5.  **Find the Minimum:**
    Our minimization problem is now:

    $$
    \min_{z' \in N(A)} \left( \|x_r\|_2^2 + \|z'\|_2^2 \right)
    $$

    The smallest possible value for $\|z'\|_2^2$ is 0, which is achieved when $z' = 0$.

6.  The unique minimum-norm solution $x^*$ occurs when $z' = 0$. This gives:

    $$
    x^* = x_r + 0 = x_r
    $$

    By definition, $x_r$ is in the row space, $x_r \in \text{span}(A^T)$. Since the row space is the orthogonal complement of the null space ($\text{span}(A^T) = N(A)^\perp$), this proves that the minimum-norm solution $x^*$ must be **orthogonal to the null space $N(A)$**.
````

## Derivation using SVD

Let's now use the **thin SVD** of $A$, which is a factorization $A = U \Sigma V^T$ where:

* $A$ is $m \times n$ with rank$(A) = r < n$.
* $U$ is an $m \times r$ matrix with orthonormal columns. Its columns form a basis for $\text{span}(A)$.
* $\Sigma$ is an $r \times r$ diagonal matrix containing the $r$ non-zero singular values, $\sigma_1 \ge \dots \ge \sigma_r > 0$. It is invertible.
* $V$ is an $n \times r$ matrix with orthonormal columns. Its columns form a basis for the **row space**, $\text{span}(A^T)$.

We will now apply our two conditions.

(ls-conditions)=
### 1. The Least-Squares Condition

The LS solution must make the residual $r = Ax - b$ orthogonal to the column space, $\text{span}(A)$.

* Since $\text{span}(A) = \text{span}(U)$, this orthogonality condition is:

    $$
    U^T (Ax - b) = 0
    $$

* Substitute $A = U \Sigma V^T$:

    $$
    U^T ( (U \Sigma V^T)x - b ) = 0
    $$

* Distribute $U^T$. Since $U$ has orthonormal columns, $U^T U = I_r$ (the $r \times r$ identity):

    $$
    \begin{aligned}
        (U^T U) \Sigma V^T x - U^T b &= 0 \\
        \Sigma V^T x &= U^T b
    \end{aligned}
    $$

    This equation defines all LS solutions. As $V^T$ is a "wide" $r \times n$ matrix, this system is underdetermined for $x$ and has infinitely many solutions.

### 2. The Minimum Norm Condition

Now we apply our second condition: the solution $x$ must be orthogonal to the null space, $x \perp N(A)$.

* From the SVD, we know the **null space $N(A)$ is the orthogonal complement of the row space $\text{span}(A^T)$**.
* Since $\text{span}(A^T) = \text{span}(V)$, the condition $x \perp N(A)$ is *equivalent* to $x \in \text{span}(V)$.
* This means $x$ must be a linear combination of the columns of $V$. We can express this by defining a new, unique vector $y \in \mathbb{R}^r$ such that:

    $$
    x = V y
    $$

### 3. Combining and Solving

We now substitute our "minimum norm" form $x = Vy$ into our "least-squares" equation from [Step 1](#ls-conditions):

$$
\Sigma V^T (Vy) = U^T b
$$

* Since $V$ has orthonormal columns, $V^T V = I_r$. We can simplify the equation:

    $$
    \begin{aligned}
        \Sigma (V^T V) y &= U^T b \\
        \Sigma y &= U^T b
    \end{aligned}
    $$

* This is now a simple $r \times r$ diagonal system. Since $\Sigma$ is invertible (it only contains non-zero singular values), we can solve for the unique $y$:

    $$
    y = \Sigma^{-1} U^T b
    $$

* Finally, we substitute $y$ back into $x = Vy$ to get our unique, minimum-norm LS solution:

    $$
    x = V y = V \Sigma^{-1} U^T b
    $$

This solution 

$$
x = \sum_{i=1}^r \frac{u_i^T b}{\sigma_i} \; v_i
$$

is the most robust and general solution to the linear least-squares problem.

:::{tip}
**The Moore-Penrose Pseudoinverse**

The matrix we constructed, $A^\dagger = V \Sigma^{-1} U^T$, is known as the **Moore-Penrose Pseudoinverse**.

The minimum-norm least-squares solution to $Ax=b$ is *always* given by:

$$
x = A^\dagger b
$$

The SVD provides the most general and numerically stable way to compute $A^\dagger$. It handles all cases:

* $A$ is invertible: $A^\dagger = A^{-1}$.
* $A$ is overdetermined (full rank): $A^\dagger = (A^T A)^{-1} A^T$.
* $A$ is rank-deficient (our case): $A^\dagger$ gives the unique minimum-norm solution.
:::

### Equivalence with the Normal Equations Solution when Full Rank

We will now prove that the minimum-norm least-squares solution using the SVD given by $V \Sigma^{-1} U^T b$ is the same as $(A^T A)^{-1} A^T b$ when $A$ has full column rank.

This proof is valid under the assumption that $A$ has **full column rank** (rank$(A) = n$), which is the same assumption required for $(A^T A)^{-1}$ to exist in the first place.

We will use the **thin SVD**, where $A = U \Sigma V^T$:

* $A$ is $m \times n$.
* $U$ is $m \times n$ with orthonormal columns, so $U^T U = I_n$.
* $\Sigma$ is $n \times n$, diagonal, and invertible (since rank$=n$, all $n$ singular values are non-zero).
* $V$ is $n \times n$ and orthogonal, so $V^T V = V V^T = I_n$ and $V^{-1} = V^T$.


````{prf:proof}
We will start with the normal equations form $A^\dagger = (A^T A)^{-1} A^T$ and substitute $A = U \Sigma V^T$ into it.

**1. Substitute $A$ into $A^T A$:**
First, let's find the expression for $A^T A$:

$$
A^T = (U \Sigma V^T)^T = V \Sigma^T U^T
$$

Since $\Sigma$ is a square diagonal matrix, $\Sigma^T = \Sigma$. So, $A^T = V \Sigma U^T$.
Now, multiply $A^T$ and $A$:

$$
A^T A = (V \Sigma U^T) (U \Sigma V^T)
$$

Group the middle terms: $A^T A = V \Sigma (U^T U) \Sigma V^T$.
Since $U^T U = I_n$ (the $n \times n$ identity), this simplifies to:

$$
A^T A = V \Sigma (I_n) \Sigma V^T = V \Sigma^2 V^T
$$

**2. Invert $A^T A$:**
Now we find the inverse of the expression from Step 1:

$$
(A^T A)^{-1} = (V \Sigma^2 V^T)^{-1}
$$

Using the inverse rule $(XYZ)^{-1} = Z^{-1} Y^{-1} X^{-1}$:

$$
(A^T A)^{-1} = (V^T)^{-1} (\Sigma^2)^{-1} V^{-1}
$$

* Since $V$ is orthogonal, $V^{-1} = V^T$ and $(V^T)^{-1} = V$.
* Since $\Sigma^2$ is diagonal, $(\Sigma^2)^{-1} = \Sigma^{-2}$.
Substituting these in, we get:

$$
(A^T A)^{-1} = V \Sigma^{-2} V^T
$$

**3. Put It All Together:**
Now we go back to the original formula and substitute our results from Step 1 and Step 2:

$$
\begin{aligned}
    A^\dagger &= (A^T A)^{-1} A^T \\
             &= (V \Sigma^{-2} V^T) (V \Sigma U^T)
\end{aligned}
$$

Group the middle terms: $A^\dagger = V \Sigma^{-2} (V^T V) \Sigma U^T$.
Since $V^T V = I_n$:

$$
A^\dagger = V \Sigma^{-2} (I_n) \Sigma U^T = V (\Sigma^{-2} \Sigma) U^T
$$

Finally, simplify the diagonal matrices: $\Sigma^{-2} \Sigma = \Sigma^{-1}$.

$$
A^\dagger = V \Sigma^{-1} U^T
$$

This completes the proof, showing the two expressions are identical.
````