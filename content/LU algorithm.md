LU factorization algorithm:

- Loop over $k$: 1 to $n$
	- $u_{k,} = a_{k,}$; row $k$ of $U$
	- $l_{,k} = a_{,k} / a_{kk}$; column $k$ of $L$
	- $A \leftarrow A - l_{,k} * u_{k,}$. In this step, we modify matrix $A$. The resulting modified matrix $A$ is called the Schur complement.

- The computational cost is $O(n^3)$.
- This algorithm assumes that $a_{kk} \neq 0$.
- $a_{kk}$ is called a **pivot** in the LU factorization. ^pivot
- This is a non-trivial assumption that can break down in practice.
- We will see later how to deal with this situation.

[[Triangular factorization]], [[Solving linear systems using LU]]