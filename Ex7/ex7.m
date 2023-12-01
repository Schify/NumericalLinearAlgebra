


function out = kernel_eval(kern, is, js, xs)
    out = zeros(length(is), lenght(js));
    for j=1:js
        for i = 1:is
        out(i, j) = kern(xs(i), xs(j));
        end
    end
end