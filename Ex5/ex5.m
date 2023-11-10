clear all
close all
n = 20;
is = (1:n)';
x = cos((2*is-1)/(2*n)*pi);
w = 1./sqrt(1-x.^2);
[D,Q] = inv_eig_prob(x, w);
figure
subplot(1,2,1)
imagesc(abs(D))
title("D")
subplot(1,2,2)
imagesc(abs(Q))
title("Q")
figure
subplot(1,2,1)
spy(D)
title("D")
subplot(1,2,2)
spy(Q)
title("Q")
function [D,Q]=inv_eig_prob(x, w)
    X = diag(x);
    D = [w X];
    N = size(X,1);
    Q = eye(N,N);
    for i=2:N
        for j = 1:(i-1)
            G = eye(N,N);
            G([j,i], [j,i])  = givens(D(i,j), D(j,j));
            Q = G*Q;
            D = G'*D*[1 zeros(1,N); zeros(N,1) G];
        end
    end
end

 function [c,s] = GivensRotation(a,b)
 %https://stackoverflow.com/questions/13438073/qr-decomposition-algorithm-using-givens-rotations
        if b == 0
            c = 1;
            s = 0;
        else
            if abs(b) > abs(a)
                r = -a / b;
                s = 1 / sqrt(1 + r^2);
                c = s*r;
            else
                r = -b / a;
                c = 1 / sqrt(1 + r^2);
                s = c*r;
            end
        end
    end