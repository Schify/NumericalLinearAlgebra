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
k= -1; 
[Q_k,T_k,r,err_ind,w_k_inf] = Lanczos2(A,k,randn(size(A,1),1));%Lanczos2
%plot_W(Q_k'*Q_k - eye(size(Q_k,2)), w_k_inf, 1900, 100, 800,250)
k = size(T_k, 1);

alpha = diag(T_k);beta=diag(T_k, 1);
eig_indicies = (1:5)'; % the smallest eigenvalues' indicies
eig_est = zeros(size(eig_indicies));
eig_vec_est = zeros(size(U,1),length(eig_est));
for j = 1:length(eig_indicies)
    eig_est(j) = eigval_k(eig_indicies(j), alpha, beta); 
    w = null(A-eig_est(j)*eye(size(A)));
    eig_vec_est(:, j)=U*S_half_inv*w;
    eig_vec_est(:, j)=eig_vec_est(:, j)./max(abs(eig_vec_est(:, j)));
end

%% Plots
figure(5)
plot_afm_mesh(nodes,elements); %plot the AFM tip
title("AFM tip mesh")
%exportgraphics(gcf,'../plots/afm_mesh_distort.pdf','ContentType','vector')

for k = 1:5
    dU = zeros(size(nodes,1)*size(nodes,2), 1);
    dU(actualDofs) = eig_vec_est(:, k);
    dU = reshape(dU, [], 3);
    displ_nodes = nodes+dU*20;
    
    figure(12+k)
    axis vis3d
    title(sprintf('Mode %i $(\\omega_%i=%g)$', k, k, sqrt(eig_est(k))), 'Interpreter', 'latex');
    plot_afm_surf(nodes, elements)
    view(3)
    hold on 
    plot_afm_frame(displ_nodes, elements)
    

    hold off
    exportgraphics(gcf,sprintf('../plots/afm_mesh_mode%i.pdf', k) ,'ContentType','vector')

end


