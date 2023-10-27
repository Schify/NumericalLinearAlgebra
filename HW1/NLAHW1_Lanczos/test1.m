addpath("../..")
close all

A=mmread("Test.mtx");

[Q_k,T_k,r,err_ind, w_k_inf]=Lanczos1(A, 25, rand(size(A,1), 1), max(max(A)));

W_k = Q_k' * Q_k;

figure("Name","Wk")
ax=imagesc(log10(abs(W_k)));
colormap("jet")
title("$W_k$ : $\log_{10} {\left| {\left[W_k\right]}_{(i,j)} \right|}$", Interpreter="latex")
xlabel("$i$", Interpreter="latex")
ylabel("$j$", Interpreter="latex")
colorbar

figure(2)
title("$W_{k, \infty}$", Interpreter="latex")
semilogy(w_k_inf)