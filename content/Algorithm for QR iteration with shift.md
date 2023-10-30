- Let's see in more detail how the [[QR iteration with shift]] works. 
- We need to use a different shift at each step. 
- How do we do this practically?

Pseudo-algorithm:

```julia
Tk = A
while not_converged:
    mu = Tk[n,n]
    Uk, Rk = qr(Tk - mu * I)
    Tk = Rk * Uk + mu * I
```

We can check that this is a valid similarity transformation using unitary matrices. Here is the sequence of steps to check that the iteration is correct and that the transformation is a unitarily similar transformation:

- $U_{k+1} R_{k+1} = T_k - \mu I$
- $R_{k+1} = U_{k+1}^H T_k - \mu U_{k+1}^H$
- $T_{k+1} = R_{k+1} U_{k+1} + \mu I$
- $T_{k+1} = (U_{k+1}^H T_k - \mu U_{k+1}^H) U_{k+1} + \mu I$
- $T_{k+1} = U_{k+1}^H \, T_k \, U_{k+1}$

So all eigenvalues are preserved, and we are converging to the $T$ matrix from the [[Schur decomposition]].