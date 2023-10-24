N = 5

A = rand(N,N); A = (A'+A)/2;
b = rand(N,1);
[Q, T] = lanczos(A, b)


function [Q, T] = lanczos(A, b)
    k = size(A, 1);
    Q = zeros(k,k+1); % Q now starts with column 0
    T = zeros(k,k);
    Q(:,2) = b / norm(b);
    alpha = zeros(k,1); beta = zeros(k,1); % need a beta0 too!!
    for j = 2:(k+1)
        z = A * Q(:,j);
        alpha(j-1) = Q(:,j)'*z;T(j-1,j-1) = alpha(j-1);
        z = z - alpha(j-1) * Q(:,j) - beta(j-1) * Q(:,j-1);
        beta(j) = norm(z);
        if beta(j) < 100*eps
            break
        end
        if j<k+1
            Q(:,j+1) = z / beta(j);

            T(j-1, j) = beta(j);
            T(j, j-1) = beta(j);
        end 
        
    end

    Q = Q(:, 2:end);

end