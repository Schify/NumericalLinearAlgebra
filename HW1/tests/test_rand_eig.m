clear all
addpath("../..")
addpath("../../matlab2tikz/src/")
addpath("../NLAHW1_Lanczos/")
addpath("../NLAHW1_AFM/")
close all
save_figs = true
test_on_rand_tridig(25, 0,save_figs)
test_on_rand_tridig(800, 0,save_figs)
test_on_rand_tridig(25, 1,save_figs)
test_on_rand_tridig(800, 1,save_figs)

function test_on_rand_tridig(k, stable,save_figs)
%% test on random matrix
    rng(1)
    % k =900;
    alpha = randn(k,1);
    beta = randn(k-1,1);
    T_k = diag(alpha)+diag(beta,1)+diag(beta,-1);
    mat_eig = eig(T_k);
    eigk = zeros(k, 1);
    % stable = true;
    for j = 1:k
        eigk(j) = eigval_k(j, alpha, beta, stable);
    end
    
    
    
    x = linspace(min(mat_eig)*1.1,max(mat_eig)*1.1,3000)';
    y = num_eigvals_smaller(x,alpha, beta, stable);
    figure
    hold on
    plot(x,y)
    
    plot(eigk, (1:k)-0.5, "x")
    %bar(mat_eig,k*ones(size(mat_eig)),0.05,"grouped", "red", "FaceAlpha",0.2);
    xline(mat_eig, "g--");
    legend(["$s(x)$", "$\mu_i$", "$\texttt{eig(A)}$"], Interpreter="latex", AutoUpdate='off', Location="southeast")
    hold off
    if  save_figs
        fprintf('eig_rand_test_k%i_s%i.tex\n',k, stable)
        cleanfigure;
%         matlab2tikz(sprintf('../plots/eig_rand_test_k%i_s%i.tex',k, stable),...
%                         'relativeDataPath','../plots',...
%                         'showInfo', true)

        set(gcf,'position',[100,100,500,400])
        exportgraphics(gcf,sprintf('../plots/eig_rand_test_k%i_s%i.pdf',k, stable),'ContentType','vector')
    end
end
