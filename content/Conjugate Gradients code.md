The [[Conjugate Gradients Version 1|CG algorithm]] is remarkably simple [[Conjugate Gradients algorithm|in the end.]] Here is the Julia code for it:

```julia
n = size(b, 1)
x = zeros(n)
r = copy(b)
p = copy(r)
rho = dot(r,r)
while sqrt(rho) > tol
    q = A*p
    mu = rho / dot(p,q) # rho = ||r_{k-1}||^2
    x += mu * p
    r -= mu * q
    rho_old, rho = rho, dot(r,r)
    p = r + (rho / rho_old) * p
end
```