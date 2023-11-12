function num_eig = num_eigvals_smaller(x, alpha, beta)
% tells how many eigenvalues the tridiagonal matrix T has less than or
% equal to x (has a size of k)
% the evaluated char poly sequence has a size of  (k, n+1)
    pn = char_polys_eval(x, alpha, beta);
    num_eig = count_flips(pn);
    %the second term is to prevent the double counting of the zeros

end

function num_flips=count_flips(pn)
% function to count the number of flips occuring in a pn sequence matrix
%respect the fact that a jth 0 inherits the sign with p_{j-1}
    n = size(pn, 2)-1;
    num_flips = sum(sign(pn(:,2:n+1)).*sign(pn(:,1:n))<0,2);%some zeros are counted double for no good reason
    num_flips = num_flips-sum(abs(sign(pn(:,2:n+1)).*sign(pn(:,1:n)))<eps,2)/2;
    %Handling those pesky zeros
%     last_sign = 1;
%     for k = 1:size(pn,1)
%         for j = 1:n
%            if abs(pn(k, j))<eps
%                last_sign = -last_sign;
%                 if sign(pn(k,j+1)) == last_sign
%                     num_flips(k) = num_flips(k)+1;
%                 else
%                     num_flips(k) = num_flips(k)+2;
%                 end
%            else
%                 last_sign = sign(pn(k,j));
%            end
%         end
%     end
end