- nnz: number of non-zero entries in matrix $A$. 
- nnz is $O(n)$ for [[Motivation of iterative methods for eigenvalue computation|sparse matrices.]]

### Arnoldi

- Matrix-vector products: $O(k \; {\rm nnz})$.
- $h_{ik} = \boldsymbol q_i^T A \boldsymbol q_k,$ $1 \le i \le k.$ $O(kn)$ per iteration.
- Vector operations: $O(k^2 n)$ total cost.

### Lanczos

- Matrix-vector products: $O(k \; {\rm nnz})$.
- Vector operations: $O(k n)$ total cost.
- Lanczos is much more efficient!

[[Algorithm for the Arnoldi process]], [[Algorithm for the Lanczos process]]