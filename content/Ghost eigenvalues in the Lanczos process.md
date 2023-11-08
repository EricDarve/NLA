- Because of roundoff errors, some eigenvalues can appear more than once. 
- This is caused by the short recurrence in [[Lanczos process|Lanczos]], as opposed to the full orthogonalization process in [[Key idea of iterative methods for eigenvalue computation|Arnoldi.]]
- This does create difficulties in practice.
- However, several modifications to [[Algorithm for the Lanczos process|Lanczos]] can be made to address this problem.

![[2022-11-05-18-37-14.png]]
![[2022-11-05-18-37-26.png]]
[[Algorithm for the Lanczos process]], [[Algorithm for the Arnoldi process]], [[Convergence of Lanczos eigenvalues for symmetric matrices]], [[Convergence of Lanczos inner eigenvalues]]