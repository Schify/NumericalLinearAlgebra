clear all
addpath("../..")
addpath("../../matlab2tikz/src/")
addpath("../NLAHW1_Lanczos/")
close all

A=mmread("Test.mtx");

%% Lanczos 3
kmax = 100;
[Q_k,T_k,r,err_ind,w_k_inf_quasi]=Lanczos3(A, kmax, rand(size(A,1), 1), max(max(A)));
kmax = size(Q_k, 2);
W_k = Q_k' * Q_k - eye(kmax, kmax);
w_k_inf = ones(1,kmax);
for k=1:kmax
    w_k_inf(k) = max(abs(W_k(1:k,1:k)),[],"all");
end


x0 = 1000; y0 = 300; width = 800; height = 250;
figure(1)

plot_W(W_k, w_k_inf_quasi, x0, y0, width, height)

cleanfigure;
matlab2tikz('../plots/lanczos3_W.tex','relativeDataPath','../plots',...
                'showInfo', false);%/lanczos0_W.tex

figure(6)
e = ones(size(T_k, 1),1);
imagesc(log10(abs(A*Q_k-Q_k*T_k)));%-r*ones(size(T_k, 1),1)' 
colorbar

real_eig_values = readmatrix("TestEig.txt");
figure(71)
plot_ritz_real(T_k, real_eig_values)
% cleanfigure;
% matlab2tikz('../plots/lanczos3_ritz.tex','relativeDataPath','../plots',...
%                 'showInfo', false)%, ...
        %'floatFormat','%.6s');%/lanczos0_W.tex
saveas(gcf,'../plots/lanczos3_ritz.pdf')

