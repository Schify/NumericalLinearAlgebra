function Q = alg2(A,l)
%Randomized Range Finder
% Given an m × n matrix A and an integer l, this scheme computes an m × l
% orthonormal matrix Q whose range approximates the range of A
[m,n] = size(A);
Omega = randn(n, l);
Y = A*Omega;
[Q, ~] = qr(Y);
end

