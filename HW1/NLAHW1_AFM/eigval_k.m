function eig_est = eigval_k(k,alpha,beta)
% finds the kth eigenvalue of a symmetric tridiagonal matrix
n = size(alpha,1);

%using the Gershgorin circle theorem to find the range of possible
%eigenvalues
%the radia
Ri = abs([0;beta(1:n-1)])+abs([beta(1:n-1);0]);%being careful at the corners!
min_eig_range = min(alpha-Ri); max_eig_range = max(alpha+Ri);

%finding the eigenvalue through interval slicing
% s(a) -> k-1 and s(b)->k
a = min_eig_range;
b = max_eig_range;
%evals = num_eigvals_smaller([a;b], alpha, beta);
%a_eval = evals(1);b_eval = evals(2);
prec = 100*eps;
while b-a > abs(b)*prec
    c = (a+b)/2;
    c_eval = num_eigvals_smaller(c, alpha, beta);
    if c_eval<k
        a = c;
    else
        b = c;
    end
end

eig_est = (a+b)/2;

end

