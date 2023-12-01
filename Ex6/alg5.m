function [U, Sigma, V] = alg5(A, Q)
%Direct SVD
% Given matrices A and Q such that (5.1) holds, this procedure computes an
% approximate factorization A ≈ U ΣV ∗ , where U and V are orthonormal, and Σ
% is a nonnegative diagonal matrix.
B = Q'*A;
[U_tilde, Sigma, V] = svd(B);
U = Q*U_tilde;
end

