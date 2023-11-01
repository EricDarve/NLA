```julia
H = zeros(kmax, kmax)
r = copy(q[:,1])
q[:,1] = r / norm(r)
for k = 1:kmax
    if k > 1
        q[:,k] = r / H[k,k-1]
    end

    r = A * q[:,k] # Multiply by A
    for i=1:k
        # Make vector orthogonal
        H[i,k] = dot(q[:,i], r)
        r -= H[i,k] * q[:,i]
    end

    if k<kmax
        H[k+1,k] = norm(r)
    end
end	
```

- This is a Gram-Schmidt like process.