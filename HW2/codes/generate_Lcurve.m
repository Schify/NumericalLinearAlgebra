function [resids, x_norms] = generate_Lcurve(solve_func,A,b, lambdas)
%generates the data needed for the L curves
    x_norms = nan * ones(size(lambdas));
    resids = nan * ones(size(lambdas));
    for i = length(lambdas) %lambdas should be in the correct orientation!!!
        x = solve_func(A, b, lambdas(i));
        x_norms(i) = norm(x);
        resids(i) = norm(A*x-b);
    end
end

