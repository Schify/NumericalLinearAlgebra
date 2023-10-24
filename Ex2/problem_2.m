function [x, res]=problem_2(A, b, alpha, k)
    [Q_tilde, H_tilde] = arnoldi(A, b, k);
    y = (eye(k,k)-alpha*H_tilde(1:k, :))\([norm(b); zeros(k-1,1)]);
    res = norm((1-alpha*H_tilde(k+1, k))*y(k,1));
    x = Q_tilde(:, 1:k)*y;
end

