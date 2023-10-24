function xsol = jacobi_method(A, b, x_init, k)
    xsol = x_init;
    M = diag(diag(A));
    N = M-A;
    xsol = iterative_method(M, N, xsol, b, k);
end