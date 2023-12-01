clear all 
close all 

load("matrices.mat")
M = M3;
tol = 1e-5;
[X, Y]  = aca(M, tol);
[U, S, V] = svd(M, "econ");
rank_of_M =  sum(diag(S)>tol);
low_rank_SVD = U(:,1:rank_of_M)*S(1:rank_of_M, 1:rank_of_M)*V(:,1:rank_of_M)';

imagesc(log(abs(M-X*Y')/max(max(abs(M))))); colorbar;

function [ik, jk]=random_pivot(B, X, Y, I, J)
    ik = randi(size(B,1)); jk = randi(size(B,2));
end