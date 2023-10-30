Let's first perform a [[QR factorization]] step using [[QR using Givens transformations|Givens]] transformations.

![[2022-10-22-16-35-11.png]]

Let's do this step-by-step.

**Step 1:** We first perform a series of small [[QR using Givens transformations|Givens transformations]] applied to the left.

![[2022-10-22-16-35-33.png]]

**Step 2:** Then these Givens transformations are applied to the right.

![[2022-10-22-16-35-52.png]]

We obtain again a matrix in upper Hessenberg form.

- Cost of QR: $O(n^2)$
- Cost of RQ: $O(n^2)$
- Total cost of one QR iteration: $O(n^2)$

This is much less than the $O(n^3)$ for the original [[QR iteration]].

### QR iteration with upper Hessenberg matrix

- The QR iteration algorithm conserves the upper Hessenberg form. 
- The time cost per iteration in upper Hessenberg form is $O(n^2)$. 
- The total time cost is $O(n^3)$.

### Algorithm:

```julia
function givens_QR_iteration_s!(H)
    n = size(H,1)
    G = zeros(2,n-1)
    for k=1:n-1
        c, s = givens(H[k,k], H[k+1,k])
        G[:,k] = [c; s]
        # Multiply by Givens rotation to the left
        apply_left_givens!(H, k, c, s)
    end
    for k=1:n-1
        # Multiply by Givens rotation to the right
        apply_right_givens!(H, k, G[1,k], G[2,k])
    end
end
```

### Summary of one step of the QR iteration with the upper Hessenberg form

![[2022-10-22-16-45-37.png]]