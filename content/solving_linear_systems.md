How to solve $Ax = b$? One of the most important computational tasks in NLA.

Notations | Description
--- | ---
$A$ | Matrix
$n$ | Size of the matrix
$A[1: k, 1: k]$ | Top left $k \times k$ block of $A$
$a_{ij}$ | entry $(i,j)$
$a_{,j}$ | column $j$
$a_{i,}$ | row $i$
det | determinant

- This section covers the use of the LU factorization to solve linear systems like $Ax = b$. 
- This method is fast and nearly optimal in terms of floating point operations. 
- However, it suffers from stability and accuracy issues in some cases.