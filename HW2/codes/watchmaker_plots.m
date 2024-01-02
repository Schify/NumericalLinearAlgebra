close all;
clear
addpath("../../../tools/regu/")
addpath("../../matlab2tikz/src/")
load("../Temp.mat");

cond(K)
n = size(K, 2);
m = size(K, 1);
Tmax = 10;
ts = linspace(0,Tmax,n)';


ps = Problem.empty;
%% Tikhinov
%%% T1
j = 1;
ps(j) = Problem(K, g, 'Tikh', "Tikhonov simple", [], ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j) = gen_data(ps(j));

%%% T2
j = j+1
L = get_l(n, 2);
ps(j) = Problem(K, g, 'Tikh', "Tikhonov $D_2$", L, ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j).Lname = "D_2";
ps(j) = gen_data(ps(j));

%%% T3
j = j+1
from_ind = 1+floor(5/10*n);
L = eye(n,n); L(from_ind:end, from_ind:end) = L(from_ind:end, from_ind:end)*100;
ps(j) = Problem(K, g, 'Tikh', "Tikhonov forced zero", L, ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j).Lname = "I_{t>5s}";
ps(j) = gen_data(ps(j));

%%% T4
j = j+1
from_ind = 1+floor(5/10*n);
L = eye(n,n); L(from_ind:end, from_ind:end) = L(from_ind:end, from_ind:end)*100;
L = [L; 300*get_l(n, 2)];
ps(j) = Problem(K, g, 'Tikh', "Tikhonov combination", L, ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j).Lname = "300D_2,I_{t>5s}";
ps(j) = gen_data(ps(j));



j_tik = j;

j_beg = 1;
j_end = j_tik;
code_name = 'tikh';
% Solution plots
fig = figure(Position=[100, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    plot(ts, f, "DisplayName","$f$ (exact solution)")
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_sols_plot.tex', code_name))

% Error plots
fig = figure(Position=[1200, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax, f)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_error_plot.tex', code_name))

% Curve plots
fig = figure(Position=[100, 700, 1000, 400]);
for i = j_beg:j_end
        ax = subplot(1,2,1);
        hold on
        plot_l_curve(ps(i), fig, ax)
        hold off 

        ax = subplot(1,2,2);
        hold on
        plot_gcv_curve(ps(i), fig, ax)
        hold off     
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_curves_plot.tex', code_name))

%% T(G)SVD
%%% TSVD
j = j+1
ps(j) = Problem(K, g, 'tsvd', "TSVD", [], ...
    (1:100));%, U, s, V);
ps(j) = gen_data(ps(j));

%%% TGSVD
j = j+1
from_ind = 1+floor(5/10*n);
L = eye(n,n); L(from_ind:end, from_ind:end) = L(from_ind:end, from_ind:end)*100;
L = [L; 300*get_l(n, 2)];
ps(j) = Problem(K, g, 'tsvd', "TGSVD", L, ...
    (1:100));%, U, s, V);
ps(j) = gen_data(ps(j));
ps(j).Lname = "300D_2,I_{t>5s}";
ps(j) = gen_data(ps(j));

j_beg = j_end+1;
j_end = j;
code_name = 'tsvd';
% Solution plots
fig = figure(Position=[100, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    plot(ts, f, "DisplayName","$f$ (exact solution)")
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_sols_plot.tex', code_name))

% Error plots
fig = figure(Position=[1200, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax, f)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_error_plot.tex', code_name))

% Curve plots
fig = figure(Position=[100, 700, 1000, 400]);
for i = j_beg:j_end
        ax = subplot(1,2,1);
        hold on
        plot_l_curve(ps(i), fig, ax)
        hold off 

        ax = subplot(1,2,2);
        hold on
        plot_gcv_curve(ps(i), fig, ax)
        hold off     
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_curves_plot.tex', code_name))

%% D(G)SVD
%%% DSVD
j = j+1
ps(j) = Problem(K, g, 'dsvd', "DSVD", [], ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j) = gen_data(ps(j));

%%% DGSVD
j = j+1
from_ind = 1+floor(5/10*n);
L = eye(n,n); L(from_ind:end, from_ind:end) = L(from_ind:end, from_ind:end)*100;
L = [L; 300*get_l(n, 2)];
ps(j) = Problem(K, g, 'dsvd', "DGSVD", L, ...
    logspace(-12, 2, 300));%, U, s, V);
ps(j) = gen_data(ps(j));
ps(j).Lname = "300D_2,I_{t>5s}";
ps(j) = gen_data(ps(j));

j_beg = j_end+1;
j_end = j;
code_name = 'dsvd';
% Solution plots
fig = figure(Position=[100, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    plot(ts, f, "DisplayName","$f$ (exact solution)")
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_sols_plot.tex', code_name))

% Error plots
fig = figure(Position=[1200, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax, f)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_error_plot.tex', code_name))

% Curve plots
fig = figure(Position=[100, 700, 1000, 400]);
for i = j_beg:j_end
        ax = subplot(1,2,1);
        hold on
        plot_l_curve(ps(i), fig, ax)
        hold off 

        ax = subplot(1,2,2);
        hold on
        plot_gcv_curve(ps(i), fig, ax)
        hold off     
end
matlab2tikz(sprintf('../plots/%s_curves_plot.tex', code_name))

%% CGLS
%%% CGLS
j = j+1
ps(j) = Problem(K, g, 'cgls', "CGLS", [], ...
    1:300);%, U, s, V);
ps(j) = gen_data(ps(j));


j_beg = j_end+1;
j_end = j;
code_name = 'cgls';
% Solution plots
fig = figure(Position=[100, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    plot(ts, f, "DisplayName","$f$ (exact solution)")
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_sols_plot.tex', code_name))

% Error plots
fig = figure(Position=[1200, 100, 1000, 400]);
for tind = 1:2 
    ax = subplot(1,2,tind);
    hold on
    for i = j_beg:j_end
        plot_best_sol(ps(i), ts, tind, fig, ax, f)
    end
    hold off 
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_error_plot.tex', code_name))

% Curve plots
fig = figure(Position=[100, 700, 1000, 400]);
for i = j_beg:j_end
        ax = subplot(1,2,1);
        hold on
        plot_l_curve(ps(i), fig, ax)
        hold off 

        ax = subplot(1,2,2);
        hold on
        plot_gcv_curve(ps(i), fig, ax)
        hold off     
end
%cleanfigure
matlab2tikz(sprintf('../plots/%s_curves_plot.tex', code_name))



printToLatexFile(ps, '../plots/trials.txt', f)





