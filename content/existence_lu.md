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

````{prf:theorem} Existence of LU
:label: thm:existence_lu_existence

The LU factorization **exists** if and only if

$$
{\rm rank}(A[1:k,1:k]) = {\rm rank}(A[1:n,1:k])
$$

for all $1 \le k \le n-1$.
````

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