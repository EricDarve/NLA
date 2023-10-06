We wish to write $A$ in the form $A = LU$ where $L$ is lower triangular and $U$ is upper triangular. 

[[Solving linear systems using LU]], [[Solving triangular systems]]

Algebraic form:
$$
a_{ij} = [LU]_{ij} = l_{i1} u_{1j} + l_{i2} u_{2j} + \dots + l_{in} u_{nj}
$$
Let's write the $LU$ product in [[Outer form of matrix-matrix product|outer form]]:
$$
A = LU = \sum_{k=1}^n l_{,k} \, u_{k,}
$$

- Column and row notations:
	- $l_{,k}$: column $k$
	- $u_{k,}$: row $k$

The factors can be computed iteratively based on the sparsity pattern of the factors and the sum decomposition.

- Let's start with column 1 in $A=LU$. 
- In $l_{,k} u_{k,}$, only $k=1$ contributes to column 1 of $A = LU$. 
- This is because entry 1 in row vector $u_{k,}$ is 0 for $k>1$.

Therefore, the equation for column 1 of $A=LU$ is
$$
a_{,1} = l_{,1} \, u_{11}
$$
- $a_{,1}$: column 1 of $A$
- $l_{,1}$: column 1 of $L$
- $u_{11}$: first entry in row vector $u_{1,}$

Looking at $a_{11}$ we get:
$$
a_{11} = l_{11} u_{11}
$$
There are many possible solutions when solving for $l_{11}$ and $u_{11}$. The simplest involves choosing $l_{11} = 1$ or $u_{11} = 1$. We choose as a convention: $l_{11} = 1$.

The final equations for column 1 are:
- $u_{1,} = a_{1,}$
- $l_{,1} = a_{,1} / a_{11}$, where we assume that $a_{11} \neq 0$.

This completely specifies the first column of $L$ and the first row of $U$.

After computing $l_{,1}$ and $u_{1,}$, we form:
$$
A - l_{,1} \, u_{1,} = \sum_{k=2}^n l_{,k} \, u_{k,}
$$
We can apply the same idea to column 2 with $k=2$. We get the 2nd column of $L$, $l_{,2}$, and 2nd row of $U$, $u_{2,}$. We then form:
$$
A - \sum_{k=1}^2 l_{,k} \, u_{k,} = \sum_{k=3}^n l_{,k} \, u_{k,}
$$
Following an iterative process with $k=3$, ..., $n$, we can compute all the columns of $L$ and rows of $U$.