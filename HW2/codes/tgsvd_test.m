close all;
clear
addpath("../../../tools/regu/")
addpath("../../matlab2tikz/src/")
load("../Temp.mat");

[U, S, V] = svd(K);
[Ug,sm,Xg,Vg] = cgsvd(K,eye(size(K)));
M=diag(sm(:,2));
k = size(K, 2);
P = zeros(k, k);
for i = 1:k
    P(i, k-i+1) = 1;
end    
spy(P)

figure("Position",[300 400 1000 323])
subplot(1,2,1)
imagesc(log10(P*abs((Vg*M)'*V)))
ax = gca;
ax.CLim = [-16, 0];
colorbar
title("$\log_{10} \left|\left(P V_g^TM V\right)_{i,j}\right| $", Interpreter="latex")
xlabel("$i$", Interpreter= "latex")
ylabel("$j$", Interpreter= "latex")
axis equal

subplot(1,2,2)
imagesc(log10(P*abs(Ug'*U)))
colorbar
title("$\log_{10} \left|\left(P U_g^TU\right)_{i,j}\right| $", Interpreter="latex")
xlabel("$i$", Interpreter= "latex")
ylabel("$j$", Interpreter= "latex")
axis equal
cleanfigure;
matlab2tikz('../plots/tgsvd_UV.tex','relativeDataPath','../plots',...
                 'showInfo', false)%, ... 'floatFormat','%.6s')


figure("Position",[1200 400 500 600])
subplot(3,1,1)
semilogy(abs(diag(S)-P*(sm(:,1)./sm(:, 2)))./abs(diag(S)))
ylabel("$\left|\frac{\gamma_{n-i}}{\mu_i} - 1  \right| $", Interpreter="latex")
subplot(3,1,2)
semilogy(abs(diag(S)-P*(sm(:,1)./sm(:, 2))))
ylabel("$\left|\gamma_{n-i}-\mu_i  \right| $", Interpreter="latex")

subplot(3,1,3)
semilogy(diag(S))
ylabel("$\mu_i$", Interpreter="latex")
xlabel("$i$", Interpreter= "latex")
cleanfigure;
matlab2tikz('../plots/tgsvd_sing_vals.tex','relativeDataPath','../plots',...
                 'showInfo', false)%, ... 'floatFormat','%.6s')

