function pns= char_polys_eval(x,alpha,beta)
% evaluate the sequence of characteristic polynomials for alpha and beta
%beta : 1...n-1
%alpha : 1..n 
%pj : 0...n -> n+1 columns of size x
n=size(alpha, 1);
pns = ones(size(x,1),n+1);%the first column is initialized with ones
pns(:,2) = alpha(1)-x;
for j = 3:n+1 % remember that pns = [p0; p1; p2; ... ; pn] while
    % beta=[beta1; beta2; ...; betan-1; (betan)] (the last element could be given but wont be used)
    pns(:,j) = (alpha(j-1)-x).*pns(:,j-1)...
                -beta(j-2)^2.*pns(:,j-2);
end

end

