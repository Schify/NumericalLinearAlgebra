%% test on random matrix
k =25;
alpha = randn(k,1);
beta = randn(k-1,1);
T_k = diag(alpha)+diag(beta,1)+diag(beta,-1);
mat_eig = eig(T_k);
eigk = zeros(k, 1);
stable = true;
for j = 1:k
    eigk(j) = eigval_k(j, alpha, beta, stable);
end



x = linspace(min(mat_eig)*1.1,max(mat_eig)*1.1,3000)';
y = num_eigvals_smaller(x,alpha, beta, stable);
figure
hold on
plot(x,y)

plot(eigk, (1:k)-0.5, "x")
xline(mat_eig, "r--");
legend(["$s(x)$", "$\mu_i$", "$\texttt{eig(A)}$"], Interpreter="latex")
hold off

