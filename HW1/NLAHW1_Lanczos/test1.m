addpath("../..")
addpath("../../matlab2tikz/src/")
close all

A=mmread("Test.mtx");

%% Lanczos 0

[Q_k,T_k,r,err_ind, w_k_inf]=Lanczos0(A, 100, rand(size(A,1), 1), max(max(A)));
k = size(Q_k, 2);
W_k = Q_k' * Q_k - eye(k, k);


x0 = 200; y0 = 100; width = 800; height = 250;
figure(1)
subplot(1,2,1)
ax=imagesc(log10(abs(W_k)));

colormap("jet")
title("$W_k$ : $\log_{10} {\left| {\left[W_k - I\right]}_{(i,j)}\right|}$", Interpreter="latex")
xlabel("$i$", Interpreter="latex")
ylabel("$j$", Interpreter="latex")
colorbar
set(gcf,'position',[x0,y0,width,height])
axis equal

subplot(1,2,2)
title("$W_{k, \infty}$", Interpreter="latex")
semilogy(w_k_inf)
xlabel("$k$", Interpreter="latex")
ylabel("$w_{k, \infty}$", Interpreter="latex")
grid on
cleanfigure;
matlab2tikz('../plots/lanczos0_W.tex','relativeDataPath','../plots');%/lanczos0_W.tex

figure(2)
imagesc(log10(abs(A*Q_k-Q_k*T_k)));%-r*ones(size(T_k, 1),1)' 
colorbar
%% Lanczos 1

[Q_k,T_k,r,err_ind, w_k_inf]=Lanczos1(A, 100, rand(size(A,1), 1), max(max(A)));

k = size(Q_k, 2);
W_k = Q_k' * Q_k - eye(k, k);

figure(3)
subplot(1,2,1)
ax=imagesc(log10(abs(W_k)));

colormap("jet")
title("$W_k$ : $\log_{10} {\left| {\left[W_k - I\right]}_{(i,j)}\right|}$", Interpreter="latex")
xlabel("$i$", Interpreter="latex")
ylabel("$j$", Interpreter="latex")
colorbar
set(gcf,'position',[x0,y0+height,width,height])
axis equal

subplot(1,2,2)
title("$W_{k, \infty}$", Interpreter="latex")
semilogy(w_k_inf)
xlabel("$k$", Interpreter="latex")
ylabel("$w_{k, \infty}$", Interpreter="latex")
grid on 
cleanfigure;
matlab2tikz('../plots/lanczos1_W.tex','relativeDataPath','../plots');%/lanczos1_W.tex

figure(4)
imagesc(log10(abs(A*Q_k-Q_k*T_k)));%-r*ones(size(T_k, 1),1)' 
colorbar