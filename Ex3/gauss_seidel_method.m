function xsol = gauss_seidel_method(A, b, x_init, k)
    xsol = x_init;
    M = diag(diag(A)) - tril(A);
    N = M-A;
    xsol = iterative_method(M, N, xsol, b, k);
end