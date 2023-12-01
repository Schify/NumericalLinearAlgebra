function val = kernel(x, y, alpha)
   d = length(x);
   js = (1:d)';
   js = reshape(js, size(x));
   gamma = 0.9.^(js-1)/pi.^alpha;
   vals = 1+(-1)^(alpha/2+1)*(2*pi)^alpha/factorial(alpha).*gamma.*bernoulli(alpha, x-y);
   val = prod(vals);
end