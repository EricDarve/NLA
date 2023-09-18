Product of matrices as a sum: $A=BC$. We can rewrite this product in the following form:
$$
a_{ij} = \sum_k b_{ik} c_{kj} = b_{i1} c_{1j} + b_{i2} c_{2j} + \dots
$$

- $b_{,k}$: column $k$ of $B$. 
- $c_{k,}$: row $k$ of $C$.

The expression above can be rewritten as:
$$
A = \sum_k b_{,k} \, c_{k,}
$$
This is the outer form of the product. This will be useful for several matrix factorizations, including the LU and QR factorizations.

[[Matrix-vector and matrix-matrix product]]