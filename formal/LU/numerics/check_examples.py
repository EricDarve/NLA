"""Numerical cross-check of the examples in content/lu_pivoting.md.

This is not part of the Lean proof. It reruns the routines of the notes in binary64 and
checks the facts that the Lean files prove for every precision:

* the 3 x 3 example: two interchanges, the packed factors, and x = (1, 0, 1)
  (`LU.A3_final`, `LU.A3_plu`, `LU.A3_solve`);
* A_eps at eps = 1e-18: fl(pi - fl(1/eps)) = -fl(1/eps), so pi is lost
  (`LU.Aeps_double`);
* the growth matrix: no interchanges, l_ij = -1, u_in = 2^(i-1), rho_n = 2^(n-1)
  (`LU.wilk_factors`, `LU.growth_wilk`);
* the optional example: y-hat, x-hat, the exact residual, and
  eta = delta / (||A||_2 / sqrt(2) + sqrt(1 + delta^2)) (`LU.solution_wilk`,
  `LU.residual_wilk`, `LU.eta_wilk`), together with floor(n/2) <= ||A||_2 <= ||A||_F <= n
  (`LU.norm2_wilk_bounds`).

It prints the backward errors of the table produced by the notes' code cell.
"""
import math

import numpy as np


def lu_factorization_with_row_pivoting(A):
    """The routine of the notes (verbatim, without the input checks)."""
    n = A.shape[0]
    P = np.eye(n)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(A[k:, k]))
        if p != k:
            A[[k, p], :] = A[[p, k], :]
            P[[k, p], :] = P[[p, k], :]
        if A[k, k] == 0:
            continue
        A[k + 1:, k] /= A[k, k]
        A[k + 1:, k + 1:] -= np.outer(A[k + 1:, k], A[k, k + 1:])
    return P


def forward_substitution(L, b):
    n = L.shape[0]
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - L[i, :i] @ y[:i]) / L[i, i]
    return y


def backward_substitution(U, y):
    n = U.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return x


def growth_matrix(n):
    A = np.eye(n) - np.tril(np.ones((n, n)), k=-1)
    A[:, -1] = 1.0
    return A


def check_3x3():
    A = np.array([[2., 1., 1.], [4., 3., 3.], [8., 7., 9.]])
    b = np.array([3., 7., 17.])
    packed = A.copy()
    P = lu_factorization_with_row_pivoting(packed)
    expected = np.array([[8, 7, 9], [1 / 4, -3 / 4, -5 / 4], [1 / 2, 2 / 3, -2 / 3]])
    assert np.allclose(packed, expected)
    assert np.array_equal(P, np.eye(3)[[2, 0, 1]])
    L = np.tril(packed, k=-1) + np.eye(3)
    U = np.triu(packed)
    y = forward_substitution(L, P @ b)
    assert np.allclose(y, [17, -5 / 4, -2 / 3])
    x = backward_substitution(U, y)
    assert np.allclose(x, [1, 0, 1])
    print("3 x 3 example: two interchanges, factors and x = (1, 0, 1) as in the notes")


def check_aeps():
    eps = 1e-18
    l21 = 1.0 / eps
    u22 = math.pi - l21 * 1.0
    assert u22 == -l21 and (0.0 - l21) == u22
    assert 2.0**59 <= l21 < 2.0**60
    print(f"A_eps, eps = 1e-18: fl(pi - fl(1/eps)) = {u22!r} = -fl(1/eps); pi is lost")


def check_growth(n):
    A = growth_matrix(n)
    packed = A.copy()
    P = lu_factorization_with_row_pivoting(packed)
    assert np.array_equal(P, np.eye(n))
    L = np.tril(packed, k=-1) + np.eye(n)
    U = np.triu(packed)
    assert np.array_equal(L, np.eye(n) - np.tril(np.ones((n, n)), k=-1))
    U_expected = np.eye(n)
    U_expected[:, -1] = 2.0 ** np.arange(n)
    assert np.array_equal(U, U_expected)
    # growth factor: largest active entry over all stages / max |a_ij|
    work = A.copy()
    rho = np.max(np.abs(work))
    for k in range(n - 1):
        lk = work[k + 1:, k] / work[k, k]
        work[k + 1:, k + 1:] -= np.outer(lk, work[k, k + 1:])
        rho = max(rho, np.max(np.abs(work[k + 1:, k + 1:])))
    assert rho / np.max(np.abs(A)) == 2.0 ** (n - 1)
    return L, U


def check_large_error():
    u = 2.0**-53
    print("\n  n   eta (notes' code)   delta/(||A||_2/sqrt2 + sqrt(1+delta^2))   floor(n/2) <= ||A||_2 <= ||A||_F <= n")
    for n in (10, 20, 30, 40, 50, 60, 70):
        A = growth_matrix(n)
        L, U = check_growth(n)
        delta = u * 2.0 ** (n - 3)
        b = np.zeros(n)
        b[0], b[-1] = 1.0, delta
        pb = b  # P = I
        y = np.zeros(n)
        for i in range(n):
            s = 0.0
            for j in range(i):
                s += L[i, j] * y[j]
            y[i] = pb[i] - s
        y_expected = np.array([1.0] + [2.0 ** (i - 1) for i in range(1, n)])
        assert np.array_equal(y, y_expected)
        xhat = backward_substitution(U, y)
        x_expected = np.zeros(n)
        x_expected[0] = x_expected[-1] = 0.5
        assert np.array_equal(xhat, x_expected)
        r = b - A @ xhat
        r_expected = np.zeros(n)
        r_expected[-1] = delta
        assert np.array_equal(r, r_expected)
        norm2 = np.linalg.norm(A, 2)
        frob = np.linalg.norm(A, "fro")
        assert n // 2 <= norm2 * (1 + 1e-12) and norm2 <= frob * (1 + 1e-12) and frob <= n
        eta = np.linalg.norm(r) / (norm2 * np.linalg.norm(xhat) + np.linalg.norm(b))
        formula = delta / (norm2 / math.sqrt(2) + math.sqrt(1 + delta**2))
        assert abs(eta - formula) <= 1e-14 * formula
        assert eta <= 1
        print(f" {n:2d}   {eta:.6e}        {formula:.6e}                              "
              f"{n // 2} <= {norm2:.4f} <= {frob:.4f} <= {n}")


if __name__ == "__main__":
    check_3x3()
    check_aeps()
    for n in range(2, 40):
        check_growth(n)
    print("growth matrix: no interchanges, l_ij = -1, u_in = 2^(i-1), rho_n = 2^(n-1) for n = 2..39")
    check_large_error()
    print("\nall checks passed")
