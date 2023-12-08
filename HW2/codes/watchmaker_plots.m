close all;
clear
addpath("../../../tools/regu/")
addpath("../../matlab2tikz/src/")
load("../Temp.mat");

cond(K)
n = size(K, 2);
m = size(K, 1);
Tmax = 10;

ps = Problem.empty;
%% Tikhinov
%%% T1
j = 1
ps(j) = Problem(K, g, 'Tikh', "Tikhonov simple", [], ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j) = gen_data(ps(j));

%%% T2
j = 2
L = get_l(n, 2);
ps(j) = Problem(K, g, 'Tikh', "Tikhonov advanced 1", L, ...
    logspace(-12, 2, 300));%, U, s, V);
proj(j).Lname = "D_2";
ps(j) = gen_data(ps(j));

%%% T3
j = 3
from_ind = 1+floor(5/10*n);
L = eye(n,n); L(from_ind:end, from_ind:end) = L(from_ind:end, from_ind:end)*100;
ps(j) = Problem(K, g, 'Tikh', "Tikhonov advanced 1", L, ...
    logspace(-12, 2, 300));%, U, s, V);
proj(j).Lname = "W_{t>5s}";
ps(j) = gen_data(ps(j));


printToLatexFile(ps, "../plots/trials.txt", f)





