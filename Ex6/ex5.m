rng(1);
A = rand(400,10)*rand(10)*rand(10,500);
algorithms = [{@alg2, @alg5};
              {@alg3, @alg5}];
len = size(algorithms,1);

[U_tsvd, S_tsvd, V_tsvd] = svd(A);
fprintf("i\t  timing\t  U_err\t\t  V_err\t\t  S_err\n-------------------------------\n")
for i = 1:len
    algs = algorithms(i,:);
    tic
    Q = algs{1}(A,15); % Stage A
    [U,S,V] = algs{2}(A, Q); % Stage B
    timing = toc();

    U_err = max(max(abs(U-U_tsvd)))/(size(U,1)*size(U,2));
    V_err = max(max(abs(V-V_tsvd)))/(size(V,1)*size(V,2));
    S_err = max(max(abs(S-S_tsvd)))/size(S,1);
    
    fprintf("%i\t%11.4s\t%11.4s\t%11.4s\t%11.4s\n", i, timing, U_err, V_err, S_err)
end