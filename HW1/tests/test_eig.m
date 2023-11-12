clear all
addpath("../..")
addpath("../NLAHW1_AFM/")
addpath("../NLAHW1_Lanczos/")
addpath("../../matlab2tikz/src/")
close all

A=mmread("../NLAHW1_Lanczos/Test.mtx");A_eigs=readmatrix("../NLAHW1_Lanczos/TestEig.txt");
%A = rand(7,7);A = (A+A')/2;A_eigs = eig(A);
kmax = 200;
[Q_k,T_k,r,err_ind,w_k_inf]=Lanczos2(A, kmax, rand(size(A,1), 1), max(max(A)));


x0 = 1000; y0 = 300; width = 800; height = 250;
figure
W_k = Q_k'*Q_k;
plot_W(W_k-eye(size(W_k)), w_k_inf, x0, y0, width, height)

alpha = full(diag(T_k));
beta = full(diag(T_k,1));
eigval_k(3, alpha,beta)

%% Test on big eigvals
% x = linspace(min(A_eigs)-1,max(A_eigs)+1,3000)';
% y = num_eigvals_smaller(x,alpha, beta );
% figure
% hold on
% plot(x,y)
% xline(A_eigs, "r--");
% hold off
% 
% figure
% T_eigs = eigs(T_k);

%% test on random matrix
k =900;
alpha = randn(k,1);
beta = randn(k-1,1);
T_k = diag(alpha)+diag(beta,1)+diag(beta,-1);
mat_eig = eig(T_k);
eigk = zeros(k, 1);
for j = 1:k
    eigk(j) = eigval_k(j, alpha, beta);
end


figure
x = linspace(min(mat_eig)*1.1,max(mat_eig)*1.1,3000)';
y = num_eigvals_smaller(x,alpha, beta );
figure
hold on
plot(x,y)

plot(eigk, (1:k)-0.5, "x")
xline(mat_eig, "r--");
legend(["$s(x)$", "$\mu_i$", "$\texttt{eig(A)}$"], Interpreter="latex")
hold off

