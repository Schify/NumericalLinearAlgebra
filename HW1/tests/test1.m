clear
addpath("../..")
addpath("../../matlab2tikz/src/")
addpath("../NLAHW1_Lanczos/")
close all

A=mmread("Test.mtx");

%% Lanczos 0

[Q_k,T_k,r,err_ind, w_k_inf]=Lanczos0(A, 100, rand(size(A,1), 1), max(max(A)));
k = size(Q_k, 2);
W_k = Q_k' * Q_k - eye(k, k);


x0 = 200; y0 = 100; width = 800; height = 250;
figure(1)
plot_W(W_k, w_k_inf, x0, y0, width, height)
cleanfigure;
matlab2tikz('../plots/lanczos0_W.tex','relativeDataPath','../plots');%/lanczos0_W.tex

figure(2)
e = ones(size(T_k, 1),1);
imagesc(log10(abs(A*Q_k-Q_k*T_k)));%-r*ones(size(T_k, 1),1)' 
colorbar

real_eig_values = readmatrix("TestEig.txt");
figure(69)
plot_ritz_real(T_k, real_eig_values)
cleanfigure;
matlab2tikz('../plots/lanczos0_ritz.tex','relativeDataPath','../plots')%, ...
        %'floatFormat','%.6s');%/lanczos0_W.tex
saveas(gcf,'../plots/lanczos0_ritz.pdf')

%% Lanczos 1

[Q_k,T_k,r,err_ind, w_k_inf]=Lanczos1(A, 100, rand(size(A,1), 1), max(max(A)));

k = size(Q_k, 2);
W_k = Q_k' * Q_k - eye(k, k);

figure(3)
plot_W(W_k,w_k_inf, x0, y0+height+100, width, height)
figure(4)
imagesc(log10(abs(A*Q_k-Q_k*T_k)));%-r*ones(size(T_k, 1),1)' 
colorbar

figure(70)
plot_ritz_real(T_k, real_eig_values)
cleanfigure;
matlab2tikz('../plots/lanczos0_ritz.tex','relativeDataPath','../plots')%, ...
        %'floatFormat','%.6s');%/lanczos0_W.tex
saveas(gcf,'../plots/lanczos1_ritz.pdf')
