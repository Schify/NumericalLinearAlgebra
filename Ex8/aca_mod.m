function [X, Y]=aca(B, tol, selection)
    %adaptive cross approximation
    epsilon = inf;
    k = 0;
    [t, s] = size(B);
    kmax = min(t, s);
    X = zeros(t, 0);
    Y = zeros(s, 0);
    I = zeros(1, 0);
    J = zeros(1, 0);
    
    if nargin < 3
        select_pivot = @partial_pivot;
    else
        select_pivot = selection;
    end

    while epsilon > tol
        [ik, jk]=select_pivot(B, X, Y, I, J);
        fprintf("%i\t%i\n", ik, jk)
        %ik = randi(size(B,1)); jk = randi(size(B,2));
        x = B(:, jk); y = B(ik, :)';
        if abs(x(ik)) < eps
            break
        end
        if k > 0
            x = x - X*Y(jk, :)';
            y = y - Y*X(ik, :)';
        end
        x = x ./ x(ik);% y = y / y(jk);
        X = [X, x]; Y = [Y, y];
        I = [I, ik]; J = [J, jk];
        epsilon = norm(x)*norm(y)/norm(X(:,1))/norm(Y(:,1));
        k = k+1;
        if k == kmax
            fprintf("reached kmax\n")
            break
        end
    end

end

function [ik, jk]=partial_pivot(B, X, Y, I, J)
    n = size(B, 2);
    available_j = setdiff((1:n)', J);
    jk = available_j(randi(n-length(J)));% jk = randi(size(B,2));
    Rk = (B(:, jk)-X*Y(jk,:)');
    y = Rk;%(jk,:);
    [~,ik] = max(abs(y));
end