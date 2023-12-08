close all;
clear
addpath("../../../tools/regu/")
addpath("../../matlab2tikz/src/")
load("../Temp.mat");

cond(K)
n = size(K, 2);
m = size(K, 1);

%[U, s, V] = csvd(K);
p1 = Problem(K, g, 'dsvd', "Tikhonov simple", [], ...
    1:200)%logspace(-12, 2, 300));%, U, s, V);
p1 = gen_data(p1);

%L = eye(m, n);
L = get_l(n, 2);
%[Ug, sm, Xg, Vg] = cgsvd(K, L);
p2 = Problem(K, g, 'tsvd', "Tikhonov advanced", L, ...
    1:200)%logspace(-12, 2, 300));%, Ug, sm, Vg, Xg);
p2 = gen_data(p2);

fig = figure;
ax = gca;
hold on
plot_gcv_curve(p1, fig, ax)
plot_gcv_curve(p2, fig, ax)
hold off


fig = figure;
ax = gca;
hold on
ts = linspace(0,10,n)';
plot_best_sol(p1,ts,2,  fig, ax)
plot_best_sol(p2,ts,2,  fig, ax)
hold off



