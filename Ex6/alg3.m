function Q = alg3(A,r, epsilon)
%Adaptive Randomized Range Finder
% Given an m × n matrix A, a tolerance ε, and an integer r (e.g., r = 10 ), the
% following scheme computes an orthonormal matrix Q such that (2) holds with
% probability at least 1 − min{m, n}*1e−r .
if nargin<3
    epsilon = 1e-6;
end

[m, n] = size(A);
Y = A*randn(n, r);

Q = zeros(m, 0);
j = 0;
while sqrt(max(sum(Y(:,j+1:j+r).^2,1)))>epsilon/(10/sqrt(2/pi))
    j = j + 1;
    if j>1
        Y(:, j) = (eye(m,m)-Q(:,1:j-1)*Q(:,1:j-1)')*Y(:,j);
        q = Y(:, j)/norm(Y(:, j));
        Q(:, j) = Q(:, 1:j-1)*q;
    else
        q = Y(:, j)/norm(Y(:, j));
        Q(:, j) = q;
    end
    
    Y(:, j+r) = (eye(m,m)-Q(:,1:j)*Q(:,1:j)')*A*randn(n,1);
    i_range = j+1:j+r-1; 
    Y(:,i_range) = Y(:,i_range)-q*(Y(:,i_range)'*q)';
end

