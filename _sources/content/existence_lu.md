# Existence and Uniqueness of LU Factorization

The LU factorization for a matrix $A$ can break down if the third **leading principal minor** is zero. This mathematical property directly causes a zero to appear on the diagonal (a zero **pivot**) during the elimination process, which halts the standard algorithm because it would require division by zero.

This process can be explored step-by-step.

## Tracing the Breakdown: A Zero Pivot Emerges 🕵️‍♀️

The LU factorization algorithm proceeds step-by-step, using the diagonal element at step $k$, known as the **pivot**, to eliminate all the entries below it. The algorithm fails if it encounters a zero pivot.

Consider performing elimination on the following example matrix:

$$
A = \begin{pmatrix}
1 & 6 & 1 & 0 \\
0 & 1 & 9 & 3 \\
1 & 6 & 1 & 1 \\
0 & 0 & 1 & 9
\end{pmatrix}
$$

**Step 1 ($k=1$):**
The first pivot is $a_{11} = 1$. It's non-zero, so the process can proceed. The $a_{31}$ entry needs to be eliminated. The multiplier is $l_{31} = a_{31}/a_{11} = 1/1 = 1$. The computation is: (Row 3) $\leftarrow$ (Row 3) - $1 \times$ (Row 1).

This gives the matrix $A^{(1)}$:

$$
A^{(1)} = \begin{pmatrix}
\mathbf{1} & 6 & 1 & 0 \\
0 & \mathbf{1} & 9 & 3 \\
0 & 0 & 0 & 1 \\
0 & 0 & 1 & 9
\end{pmatrix}
$$

**Step 2 ($k=2$):**
The second pivot is $a_{22}^{(1)} = 1$. It's also non-zero. The entries below it in the second column are already zero, so no elimination is needed.

**Step 3 ($k=3$):**
The third pivot is $a_{33}^{(2)} = \mathbf{0}$. To continue, this pivot would be used to eliminate the entry $a_{43}^{(2)} = 1$. This would require calculating the multiplier $l_{43} = a_{43}^{(2)}/a_{33}^{(2)} = 1/0$.

This is an undefined operation. The algorithm cannot proceed. The matrix $U$ would have a zero on its diagonal ($u_{33}=0$), and it cannot be used to complete the elimination.

## LU and Determinants

We will need the following fact in the proof of the theorem below.

One of the most elegant applications of LU factorization is computing the determinant. Given $A=LU$:

$$\det(A) = \det(L)\det(U)$$

* Since $L$ is unit lower triangular, $\det(L) = 1$.
* Since $U$ is upper triangular, its determinant is the product of its diagonal entries (the pivots).

Therefore:

$$\det(A) = \prod_{i=1}^n u_{ii}$$

This reduces the computationally expensive task of calculating a determinant to the cost of performing an LU factorization and then multiplying $n$ numbers.

## The Condition for Existence: Leading Principal Minors

The breakdown is not a coincidence; it is a direct consequence of a fundamental property of the matrix $A$. The existence of a standard LU factorization is completely determined by the **leading principal submatrices** of $A$.

The $k$-th leading principal submatrix, denoted as $A_k$, is the $k \times k$ matrix in the upper-left corner of $A$. Its determinant, $\det(A_k)$, is the $k$-th **leading principal minor**.

### **Theorem: Existence and Uniqueness of LU**

````{prf:theorem} Uniqueness and Existence of LU Factorization
:label: thm:existence_lu

An $n \times n$ matrix $A$ has a **unique** LU factorization (where $L$ has 1s on its diagonal) if and only if all its leading principal minors are non-zero.

$$\det(A_k) \neq 0 \quad \text{for } k = 1, 2, \dots, n-1.$$
````

Why is this true? The proof reveals a beautiful connection between the pivots and these minors. From the matrix partition $A_k = L_k U_k$, we can take the determinant:

$$\det(A_k) = \det(L_k) \det(U_k)$$

Since $L_k$ is unit triangular, its determinant is 1. The determinant of the upper triangular $U_k$ is the product of its diagonal entries:

$$\det(A_k) = 1 \cdot (u_{11} u_{22} \cdots u_{kk})$$

From this, any pivot $u_{kk}$ can be isolated:

$$u_{kk} = \frac{u_{11} u_{22} \cdots u_{kk}}{u_{11} u_{22} \cdots u_{k-1,k-1}} = \frac{\det(A_k)}{\det(A_{k-1})} \quad (\text{defining } \det(A_0) = 1)$$

This formula is the key! It shows that the $k$-th pivot, $u_{kk}$, is determined entirely by the ratio of two consecutive leading principal minors. The pivot will be well-defined and non-zero if and only if $\det(A_{k-1})$ and $\det(A_k)$ are both non-zero. For the entire algorithm to succeed (at least up to the last step), it is required that $\det(A_k) \neq 0$ for $k=1, \dots, n-1$.

### Applying the Theorem to the Example

Let's compute the leading principal minors of the example matrix $A$:
* **$k=1$**: $\det(A_1) = \det(1) = 1 \neq 0$.
* **$k=2$**: $\det(A_2) = \det \begin{pmatrix} 1 & 6 \\ 0 & 1 \end{pmatrix} = 1 \neq 0$.
* **$k=3$**: $\det(A_3) = \det \begin{pmatrix} 1 & 6 & 1 \\ 0 & 1 & 9 \\ 1 & 6 & 1 \end{pmatrix} = 0$. (Note: subtract Row 1 from Row 3, which gives a row of zeros, making the determinant clearly 0).

The theorem predicts the factorization will fail when computing the third pivot, since $\det(A_3)=0$. Using the formula:

$$u_{33} = \frac{\det(A_3)}{\det(A_2)} = \frac{0}{1} = 0$$

This confirms precisely what was observed during the elimination process.

## Generalizations

### What if a matrix is singular?

Notice the theorem only requires the minors to be non-zero up to $k=n-1$. What if $\det(A_n) = \det(A) = 0$? The matrix is singular. The factorization can still exist, but the final pivot, $u_{nn} = \det(A_n)/\det(A_{n-1})$, will be zero. This is perfectly fine, as no further eliminations depend on $u_{nn}$.

### Existence

Here is a more precise existence statement:

```{prf:theorem} Existence of LU
:label: thm:existence_lu_existence

The LU factorization **exists** if and only if

$$
{\rm rank}(A[1:k,1:k]) = {\rm rank}(A[1:n,1:k])
$$

for all $1 \le k \le n-1$.
```

<!--

```{prf:proof}

**Notation.** To simplify the proof, let's define two submatrices:

* $A_k = A[1:k, 1:k]$: the $k$-th leading principal submatrix of $A$.
* $B_k = A[1:n, 1:k]$: the first $k$ columns of $A$.

Note that $A_k$ consists of the first $k$ rows of $B_k$. We can write $B_k$ in block form as:

$$
B_k = \begin{bmatrix} A_k \\ A_{21} \end{bmatrix},
$$ 

where $A_{21} = A[k+1:n, 1:k].$

The theorem can thus be restated as: 

$$
\boxed{\text{$A=LU$ exists $\iff \text{rank}(A_k) = \text{rank}(B_k)$ for $k=1, \dots, n-1$.}}
$$

The condition 

$$
\text{rank}(A_k) = \text{rank}\left(\begin{bmatrix} A_k \\ A_{21} \end{bmatrix}\right)
$$

is equivalent to stating that the rows of $A_{21}$ are in the row space of $A_k$.

The proof has two parts: necessity ($\implies$) and sufficiency ($\impliedby$).

**Necessity ($\implies$)**

Assume an LU factorization $A = LU$ exists. We must show that $\text{rank}(A_k) = \text{rank}(B_k)$ for all $k=1, \dots, n-1$.

Let's partition $A, L, U$ at the $k$-th row and column:

$$
A = \begin{bmatrix} A_k & A_{12} \\ A_{21} & A_{22} \end{bmatrix} = L U = \begin{bmatrix} L_k & 0 \\ L_{21} & L_{22} \end{bmatrix} \begin{bmatrix} U_k & U_{12} \\ 0 & U_{22} \end{bmatrix}
$$

Here, $L_k = L[1:k, 1:k]$ and $U_k = U[1:k, 1:k]$. $L_k$ is unit lower triangular.

From this block multiplication, we can extract two equations:

$$
\begin{aligned}
A_k &= L_k U_k \\
B_k &= A[1:n, 1:k] = \begin{bmatrix} A_k \\ A_{21} \end{bmatrix} = \begin{bmatrix} L_k U_k \\ L_{21} U_k \end{bmatrix} = \begin{bmatrix} L_k \\ L_{21} \end{bmatrix} U_k
\end{aligned}
$$

Let's analyze the ranks of these expressions.

* **From (1):** $L_k$ is a $k \times k$ unit lower triangular matrix, which means $\det(L_k) = 1$. Therefore, $L_k$ is invertible. Since we are multiplying $U_k$ by an invertible matrix, the rank is preserved:

    $$
    \text{rank}(A_k) = \text{rank}(L_k U_k) = \text{rank}(U_k)
    $$

* **From (2):** Let 

    $$
    C_k = \begin{bmatrix} L_k \\ L_{21} \end{bmatrix}.
    $$
    
    This $n \times k$ matrix is composed of the first $k$ columns of $L$. Since $L$ is unit lower triangular, its $j$-th column has a $1$ in the $j$-th position and zeros above it. This structure guarantees that the first $k$ columns of $L$ are linearly independent.
    Thus, $C_k$ has full column rank, i.e., $\text{rank}(C_k) = k$.
    Multiplying $U_k$ on the left by a matrix with full column rank preserves the rank of $U_k$:

    $$
    \text{rank}(B_k) = \text{rank}(C_k U_k) = \text{rank}(U_k)
    $$

Combining these two results, we have $\text{rank}(A_k) = \text{rank}(U_k)$ and $\text{rank}(B_k) = \text{rank}(U_k)$.

Therefore, $\text{rank}(A_k) = \text{rank}(B_k)$. This holds for all $k=1, \dots, n-1$.

**Sufficiency ($\impliedby$)**

This is a challenging part. We assume $\text{rank}(A_k) = \text{rank}(B_k)$ for $k=1, \dots, n-1$ and show $A=LU$ exists by induction on the matrix size $n$.

**Base Case (n=1):** $A = [a_{11}]$. We can choose $L=[1]$ and $U=[a_{11}]$. The rank condition for $k=1, \dots, n-1$ is empty. The statement holds.

**Inductive Step:** Assume the theorem holds for all $(n-1) \times (n-1)$ matrices. We must show it holds for an $n \times n$ matrix $A$.
We partition $A$ and seek a corresponding $L$ and $U$:

$$
A = \begin{bmatrix} a_{11} & \mathbf{w}^T \\ \mathbf{v} & C \end{bmatrix} = \begin{bmatrix} 1 & \mathbf{0}^T \\ \mathbf{l}_{21} & L' \end{bmatrix} \begin{bmatrix} u_{11} & \mathbf{u}_{12}^T \\ \mathbf{0} & U' \end{bmatrix}
$$

where $C, L', U' \in \mathbb{R}^{(n-1) \times (n-1)}$.
This requires:
1.  $u_{11} = a_{11}$
2.  $\mathbf{u}_{12}^T = \mathbf{w}^T$
3.  $\mathbf{l}_{21} u_{11} = \mathbf{v} \implies \mathbf{l}_{21} a_{11} = \mathbf{v}$
4.  $A' = L'U'$, where $A' = C - \mathbf{l}_{21} \mathbf{u}_{12}^T = C - \mathbf{l}_{21} \mathbf{w}^T$

We must find $\mathbf{l}_{21}$ and show that the resulting Schur complement $A'$ *also* satisfies the rank condition, allowing us to apply the inductive hypothesis.

* **Step 1: Find $\mathbf{l}_{21}$**
    The rank condition for $k=1$ is $\text{rank}(A_1) = \text{rank}(B_1)$, which means $\text{rank}([a_{11}]) = \text{rank}(\begin{bmatrix} a_{11} \\ \mathbf{v} \end{bmatrix})$.
    * **Case 1: $a_{11} \neq 0$.** $\text{rank}([a_{11}]) = 1$. This implies $\text{rank}(\begin{bmatrix} a_{11} \\ \mathbf{v} \end{bmatrix}) = 1$, which means $\mathbf{v}$ is a multiple of $a_{11}$. We can **uniquely** choose $\mathbf{l}_{21} = \mathbf{v} / a_{11}$.
    * **Case 2: $a_{11} = 0$.** $\text{rank}([a_{11}]) = 0$. The condition forces $\text{rank}(\begin{bmatrix} 0 \\ \mathbf{v} \end{bmatrix}) = 0$, which means $\mathbf{v} = \mathbf{0}$. The equation $\mathbf{l}_{21} \cdot 0 = \mathbf{0}$ is satisfied by **any** vector $\mathbf{l}_{21}$.

* **Step 2: Show $A'$ satisfies the rank condition**
    This is the core of the proof. We must show there *exists* an $\mathbf{l}_{21}$ (from Case 1 or 2) such that $A' = C - \mathbf{l}_{21}\mathbf{w}^T$ satisfies $\text{rank}(A'_k) = \text{rank}(B'_k)$ for $k=1, \dots, n-2$.

    * **Case 1: $a_{11} \neq 0$.**
        We chose the unique $\mathbf{l}_{21} = \mathbf{v}/a_{11}$. Let $L_1 = \begin{bmatrix} 1 & \mathbf{0}^T \\ \mathbf{l}_{21} & I \end{bmatrix}$. $L_1$ is invertible.
        Let $A^{(2)} = L_1^{-1} A = \begin{bmatrix} 1 & \mathbf{0}^T \\ -\mathbf{l}_{21} & I \end{bmatrix} \begin{bmatrix} a_{11} & \mathbf{w}^T \\ \mathbf{v} & C \end{bmatrix} = \begin{bmatrix} a_{11} & \mathbf{w}^T \\ \mathbf{0} & C - \mathbf{l}_{21}\mathbf{w}^T \end{bmatrix} = \begin{bmatrix} a_{11} & \mathbf{w}^T \\ \mathbf{0} & A' \end{bmatrix}$.
        Since $L_1^{-1}$ and its principal submatrices are invertible, they preserve the ranks of the submatrices of $A$. So $A$ satisfies the rank condition iff $A^{(2)}$ does.
        For $k \ge 1$, let $k' = k+1$. The rank condition for $A^{(2)}$ at $k'$ is:
        $\text{rank}(A^{(2)}_{k'}) = \text{rank}(B^{(2)}_{k'})$
        $\text{rank}\left(\begin{bmatrix} a_{11} & \mathbf{w}_k^T \\ \mathbf{0} & A'_k \end{bmatrix}\right) = \text{rank}\left(\begin{bmatrix} a_{11} & \mathbf{w}_k^T \\ \mathbf{0} & B'_k \end{bmatrix}\right)$
        Because $a_{11} \neq 0$, the rank of these block-upper-triangular matrices is $1 + \text{rank}(\text{Schur complement})$. By performing row operations (which preserve rank), this is equivalent to:
        $\text{rank}\left(\begin{bmatrix} a_{11} & \mathbf{0}^T \\ \mathbf{0} & A'_k \end{bmatrix}\right) = \text{rank}\left(\begin{bmatrix} a_{11} & \mathbf{0}^T \\ \mathbf{0} & B'_k \end{bmatrix}\right)$
        $1 + \text{rank}(A'_k) = 1 + \text{rank}(B'_k)$
        This implies $\text{rank}(A'_k) = \text{rank}(B'_k)$ for all $k=1, \dots, n-2$. The inductive hypothesis holds.

    * **Case 2: $a_{11} = 0$ and $\mathbf{v} = \mathbf{0}$.**
        Here $A = \begin{bmatrix} 0 & \mathbf{w}^T \\ \mathbf{0} & C \end{bmatrix}$. This is where the simple proof fails. It is **not** guaranteed that the choice $\mathbf{l}_{21} = \mathbf{0}$ will work. However, it can be proven (as by Okunev and Johnson) that the rank condition on $A$ is precisely the necessary and sufficient condition to guarantee that **there exists** a vector $\mathbf{l}_{21}$ such that the Schur complement $A' = C - \mathbf{l}_{21}\mathbf{w}^T$ *does* satisfy the rank condition.

        Let's see why this is plausible (a full proof is very technical):
        * If $\mathbf{w} = \mathbf{0}$, then $A = \begin{bmatrix} 0 & \mathbf{0}^T \\ \mathbf{0} & C \end{bmatrix}$. We can choose $\mathbf{l}_{21} = \mathbf{0}$, so $A' = C$. The rank condition on $A$ for $k \ge 2$ simplifies to $\text{rank}(C_{k-1}) = \text{rank}(D_{k-1})$, which is exactly the rank condition required for $C$. The induction holds.
        * If $\mathbf{w} \neq \mathbf{0}$, we must find a non-zero $\mathbf{l}_{21}$. The rank condition on $A$, $\text{rank}(\begin{bmatrix} \mathbf{w}_{k-1}^T \\ C_{k-1} \end{bmatrix}) = \text{rank}(\begin{bmatrix} \mathbf{w}_{k-1}^T \\ D_{k-1} \end{bmatrix})$, provides a non-trivial constraint on the relationship between $\mathbf{w}$ and $C$. This constraint is strong enough to prove that a vector $\mathbf{l}_{21}$ exists that "repairs" the rank condition for the sub-problem $A'$.

Since in all cases (Case 1 and Case 2) we can find an $\mathbf{l}_{21}$ that leads to a Schur complement $A'$ satisfying the theorem's conditions, the inductive hypothesis applies. $A' = L'U'$ exists.

We can then construct $L$ and $U$ for $A$, completing the induction.

By induction, the theorem holds for all $n$.
```

-->

## Numerical Instability 📉

Small pivots in LU factorization are numerically disastrous because they introduce massive numbers into the $L$ and $U$ factors. In finite-precision computer arithmetic, these large numbers can completely erase the original, smaller-scale information in the matrix, a phenomenon known as **catastrophic cancellation**. This leads to a computed factorization that is wildly inaccurate.

A core principle in numerical analysis is that if an algorithm fails for a specific input, it will often produce large errors for inputs that are *close* to that failure point.

We've already seen that the LU algorithm fails when a pivot is exactly zero (e.g., $a_{11}=0$). This leads us to a crucial question: what happens if a pivot is not exactly zero, but just very small?

Let's investigate with your example. The matrix

$$
A =
\begin{pmatrix}
0 & 1 \\
1 & \pi
\end{pmatrix}
$$

has no LU factorization because $a_{11}=0$. Now, consider a matrix that is nearly identical by letting the top-left entry be a very small, non-zero number $\epsilon$.

### The Small Pivot Problem in Action 🔬

Let's analyze the LU factorization of the perturbed matrix:

$$
A =
\begin{pmatrix}
    \epsilon & 1 \\
    1 & \pi
\end{pmatrix}
$$

The exact factorization is:

$$
L = \begin{pmatrix}
    1 & 0 \\ 1/\epsilon & 1
\end{pmatrix}, \quad
U = \begin{pmatrix}
    \epsilon & 1 \\ 0 & \pi - 1/\epsilon
\end{pmatrix}
$$

The first sign of trouble is immediate: a very small number, $\epsilon$, in the original matrix $A$ has created very large numbers, $1/\epsilon$, in the factors $L$ and $U$. This amplification of magnitudes is a classic red flag for numerical instability.

### Catastrophic Cancellation in Floating-Point Arithmetic

Computers cannot store real numbers with infinite precision. They use **floating-point arithmetic**, which typically keeps about 16 significant decimal digits for any number (this is the IEEE 754 double-precision standard).

Now, let's examine the term $u_{22} = \pi - 1/\epsilon$.

If $\epsilon$ is very small, say $\epsilon = 10^{-18}$, then $1/\epsilon = 10^{18}$. The computer must calculate:

$$3.1415926535... - 1,000,000,000,000,000,000$$

To perform this subtraction, the computer must align the decimal points. All the significant digits of $\pi$ are shifted so far to the right that they fall outside the 16-digit precision window. The smaller number is completely overwhelmed by the larger one.

This effect is called **swamping** or **catastrophic cancellation**. The information contained in $\pi$ is completely lost. The computer effectively calculates:

$$\pi - 1/\epsilon \approx -1/\epsilon$$

### The Corrupted Result

This single roundoff error corrupts the entire factorization. The computed upper triangular matrix, let's call it $\tilde{U}$, becomes:

$$
\tilde{U} = \begin{pmatrix}
\epsilon & 1 \\ 0 & -1/\epsilon
\end{pmatrix}
$$

If we now multiply our computed factors $\tilde{L}$ (which is just $L$) and $\tilde{U}$ back together, we don't get our original matrix $A$:

$$
\tilde{L}\tilde{U} =
\begin{pmatrix}
    1 & 0 \\ 1/\epsilon & 1
\end{pmatrix}
\begin{pmatrix}
    \epsilon & 1 \\ 0 & -1/\epsilon
\end{pmatrix}
=
\begin{pmatrix}
    \epsilon & 1 \\ 1 & 0
\end{pmatrix}
\neq
\begin{pmatrix}
    \epsilon & 1 \\
    1 & \pi
\end{pmatrix}
$$

The original $\pi$ in the matrix has vanished and been replaced by a 0. Any solution to a linear system $Ax=b$ based on this incorrect factorization will be completely wrong.