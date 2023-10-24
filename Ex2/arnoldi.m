function [Q, H] = arnoldi (A, b, k)
    if nargin<3
        k = min(size(A, 1),50);
    end
    n = size(A,1);
    Q = zeros(n,k+1); H = zeros(k+1,k);
    Q(:,1) = b / norm(b);
    for j = 1:k
        z = A*Q(:,j);
        for i = 1:j
%             H(i,j) = Q(:,i)'*z;%This seems to not work in the double
%             % orthogonolization case
            z = z - H(i,j)*Q(:,i);
        end
        for i = 1:j
            H(i,j) = Q(:,i)'*z;
            z = z - H(i,j)*Q(:,i);
        end

        if(j<k+1)%not to index out of the matrix
            H(j+1, j) = norm(z);
    
            if H(j+1, j) == 0 %Lucky breakdown
                break
            end
    
            Q(:,j+1) = z/H(j+1,j);
        end
    end

end