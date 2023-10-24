function x_next=iterative_method(M, N, x, b, k) 
    x_next = x;
    for i = 1:k
        x_next = M\(N*x+b);
    end
end