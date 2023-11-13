clear all
addpath("../..")
addpath("../../matlab2tikz/src/")
addpath("../NLAHW1_Lanczos/")
addpath("../NLAHW1_AFM/")
close all
load("../NLAHW1_AFM/HW1.mat")
%% *Converting general eigenvalue problem into standard*
% the two matrices are approximately orthogonal, so we wanna preserve that
% for later

% we perform cholesky factorization via svd
[U, S, V] = svd(full(M+M')); % $M = U S V^\intercal$, could potentially be sparse step
S_half_inv = diag(sqrt(diag(S)).^-1);
A = S_half_inv * U' * (K+K') * V *S_half_inv; % A is still symmertic,
% and now we have:
% $U S^{\frac{1}{2}} \bigl( \underbrace{S^{-\frac{1}{2}} U^\intercal K V S^{-\frac{1}{2}}}_{A} - \omega^2 I\bigr)S^{\frac{1}{2}} V^\intercal \delta \mathbf{u} = 0 $

%% *Solving the eigenvalue problem with Lanczos*
k= 300; 
[Q_k,T_k,r,err_ind,w_k_inf] = Lanczos2(A,k,randn(size(A,1),1));%Lanczos2
plot_W(Q_k'*Q_k - eye(size(Q_k,2)), w_k_inf, 1900, 100, 800,250)
k = size(T_k, 1);

alpha = diag(T_k);beta=diag(T_k, 1);
eig_indicies = (1:60)'; % the smallest eigenvalues' indicies
eig_est = zeros(size(eig_indicies)); 
for j = 1:length(eig_indicies)
    eig_est(j) = eigval_k(eig_indicies(j), alpha, beta); 
end

%%Checking against matlab eig fucntion
mat_eigs = sort(readmatrix("../NLAHW1_AFM/matlab_eigs.csv"));
mat_eigs_A = sort(readmatrix("../NLAHW1_AFM/matlab_A_eigs.csv"));

figure(25)
clf
hold on 
semilogy(mat_eigs)
semilogy(mat_eigs_A)
semilogy((1:k)./k*length(mat_eigs),eig(T_k), "x")
semilogy((eig_indicies')./k*length(mat_eigs),eig_est', "x")
set(gca, "YScale", "log")
legend(["$\texttt{eig(K,M)}$", "$\texttt{eig(A)}$","$\texttt{eig(T\_k)}$",...
    "$\texttt{Lanczos2} \, \& \,\texttt{eigval\_k}$"], Interpreter="latex", ...
    Location="southeast")
grid on
xlabel("Index of the eigenvalue (max at the full spectrum)")
ylabel("Eigenvalue")
hold off
        cleanfigure;
        matlab2tikz('../plots/afm_eigs.tex',...
                    'relativeDataPath','../plots',...
                    'showInfo', false)
