- Consider a general linear system $Ax=b$.
- Assume we have $A = LU$ where $L$ is lower triangular and $U$ is upper triangular.

We can solve $Ax = b$ in two steps:
1. $Lz = b$
2. $Ux = z$

Both steps involve [[Solving triangular systems|triangular systems]].

The computational cost of solving $Ax=b,$ assuming we have $A$ in the form $A=LU,$ is $O(n^2)$.

How can we get $L$ and $U$ starting from $A$?