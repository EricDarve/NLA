The previous methods ([[Jacobi iteration|Jacobi]], [[Gauss-Seidel iteration|Gauss-Seidel]]) lack any parameter that can be adjusted to accelerate convergence in difficult cases.

Let us revisit [[Gauss-Seidel iteration|Gauss-Seidel:]]
$$
D x^{(k+1)}_{GS} = b + L x^{(k+1)}_{GS} + U x^{(k)}
$$
Let us rewrite this as
$$
x^{(k+1)}_{GS} = x^{(k)} + \Delta x_{GS}
$$
Then SOR is defined as
$$
x^{(k+1)}_{SOR} = x^{(k)} + \omega \, \Delta x_{GS}
$$
for some parameter $\omega$.

### SOR: Successive Over-Relaxation

[[Gauss-Seidel iteration|Gauss-Seidel]] update:
$$
x^{(k+1)}_{GS} = x^{(k)} + \bigg[ D^{-1}
\big( b + L x^{(k+1)}_{GS} + U x^{(k)} \big) - x^{(k)} \bigg]
$$
Here is the increment form of this update:
$$
\Delta x_{GS} =
D^{-1} \big( b + L x^{(k+1)}_{GS} + U x^{(k)} \big) - x^{(k)}
$$
We can now define the SOR update:
$$
x^{(k+1)}_{SOR} = x^{(k)}_{SOR} + \omega \bigg[ D^{-1}
\big( b + L x^{(k+1)}_{SOR} + U x^{(k)}_{SOR} \big) - x^{(k)}_{SOR} \bigg]
$$

## SOR relaxation parameter

$\omega$ can be used for different purposes.

$$
x^{(k+1)}_{SOR} = x^{(k)}_{SOR} + \omega \, \Delta x
$$

- If the method is diverging, we can apply a smaller increment: $\omega < 1$.
- $\omega > 1$: when possible, we can try to accelerate convergence by applying a larger increment.
- This method is a little similar to the gradient descent optimization algorithm with momentum, if you are familiar with this approach.

## Code for SOR

```julia
for k = 1:maxiter
    for i = 1:n
        sigma = 0.
        for j = 1:n
            if j != i
                sigma += A[i,j] * x[j]
            end
        end
        sigma = (b[i] - sigma) / A[i,i]
        x[i] += omega * (sigma - x[i])
    end
end
```

## Ordering

As in [[Gauss-Seidel iteration|Gauss-Seidel,]] the ordering matters in this algorithm. Here are 3 simple orderings that are commonly used:

1. 1 to $n$
2. 1 to $n$ followed by $n$ to 1 in **Symmetric SOR**
3. Or, a different randomized permutation at every iteration