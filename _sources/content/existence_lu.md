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