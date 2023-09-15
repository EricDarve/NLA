How to solve $Ax = b$? One of the most important computational tasks in NLA.

## Solving triangular systems

Solve triangular system.

Draw system and solution process.

## General systems

Assume we have $A = LU$.

Show how to solve in 2 steps.

## Triangular factorization

How can we get $L$ and $U$?

Product of matrices as a sum.
$$
A = BC\\
a_{ij} = \sum_k b_{ik} c_{kj} = b_{i1} c_{1j} + b_{i2} c_{2j} + \dots\\
A = \sum_k b_{,k} \, c_{k,}
$$

## Application to LU
$$
[LU]_{ij} = l_{i1} u_{1j} + l_{i2} u_{2j} + \dots + l_{in} u_{nj}\\[1em]
LU = \sum_k l_{,k} \, u_{k,}
$$
- Column and row notations:
	- $a_{i,}$: row $i$
	- $a_{,j}$: column $j$

Compact notation:
$$
LU = \sum_k l_{,k} \; u_{k,}
$$

Explain the process of computing the factors based on the sparsity pattern of the factors and the sum decomposition.

Equations for column 1:
$$
a_{,1} = l_{,1} \; u_{11}
$$

Solution with $l_{11} = 1$ or $u_{11} = 1$. Choose: $l_{11} = 1$.

- Final equations:
	- $u_{1,} = a_{1,}$
	- $l_{,1} = a_{,1} / a_{11}$

**Full algorithm**

- Loop over $k$: 1 to $n$:
	- $u_{k,} = a_{k,}$
	- $l_{,k} = a_{,k} / a_{kk}$
	- $A \leftarrow A - l_{,k} * u_{k,}$
	- This is called the Schur complement.
