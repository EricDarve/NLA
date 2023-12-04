Here is a short Julia code implementing this [[Preconditioned Conjugate Gradients algorithm|algorithm:]]

```julia
n = size(b, 1)
x = zeros(n)
r = copy(b)
z = M*r
p = copy(z)
rho = dot(r,z)
while sqrt(rho) > tol
    q = A*p
    mu = rho / dot(p,q) # rho = ||r_{k-1}||^2
    x += mu * p
    r -= mu * q
    z = M*r
    rho_old, rho = rho, dot(r,z)
    p = z + (rho / rho_old) * p
end
```

It is very similar to the [[Conjugate Gradients code|CG code.]]