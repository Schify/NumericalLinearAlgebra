close all

err = zeros(200,1);
for N = 1:size(err, 1)
A = rand(N,N);
b = rand(N,1);
[Q, H] = arnoldi(A, b);
err(N) = norm(Q'*A*Q-H,"fro");
end

figure
loglog(err)


function [Q, H]=arnoldi(A, b)
    k = size(A, 1);
    Q = zeros(k,k); H = zeros(k,k);
    Q(:,1) = b / norm(b);
    for j = 1:k
        z = A*Q(:,j);
        for i = 1:j
%             H(i,j) = Q(:,i)'*z;%This seems to not work in the double
%             orthogonolization case
            z = z - H(i,j)*Q(:,i);
        end
        for i = 1:j
            H(i,j) = Q(:,i)'*z;
            z = z - H(i,j)*Q(:,i);
        end

        if(j<k)%not to index out of the matrix
            H(j+1, j) = norm(z);
    
            if H(j+1, j) == 0 %Lucky breakdown
                break
            end
    
            Q(:,j+1) = z/H(j+1,j);
        end
    end

end