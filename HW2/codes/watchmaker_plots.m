close all;
clear
addpath("../../../tools/regu/")
addpath("../../matlab2tikz/src/")
load("../Temp.mat");

cond(K)
n = size(K, 2);
m = size(K, 1);

[U, s, V] = csvd(K);
p1 = Problem(K, g, 'Tikh', "Tikhonov simple", [], ...
    logspace(-13, -1, 300), U, s, V);
p1 = gen_data(p1);
%loglog(p1.eta, p1.rho)
plot_l_curve(p1)

