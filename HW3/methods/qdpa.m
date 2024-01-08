function [lambda_hat, x, y] = qdpa(M, D, K, b, c, s0, tol)
%DPA Dominant Pole Algorithm
% Inputs: system : sEx = Ax+b
%                   y = C^T x
%         s0     : initial s guess (must be complex)
% Outputs: lambda_hat : the \hat{\lambda} dominant pole,
%          x, y       : eigenpair (right and left eigenvector)
k = 0;
err = inf;
sk = s0;
temp = (sk^2*M+sk*D+K);
while err>tol
    vk = temp\b;
    wk = temp'\c; 
    
    % compute new pole estimate:
    sk = sk - c'*vk/(wk'*(2*sk*M+D)*vk);
    x = vk/norm(vk);
    y = wk/norm(wk); % might not be necessary each iteration
    temp =(sk^2*M+sk*D+K);
    err = norm(temp*x);
    k=k+1;
end
lambda_hat = sk;
%x = vk; y = wk;
end

