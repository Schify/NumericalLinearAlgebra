clear 
close all 
addpath("../matlab2tikz/src/")
addpath("./methods/")
addpath("../")

N=50;
%% 1 
omega_exp_range = [0, 10];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
[A, B, E] = read_system(1);
H_og = @(s) B'*((E*s-A)\B);
bode_vals = bode_from_system(A, E, B, B, s_vals);
legend_list = ["$H(s)$"];
plot_bode(omega_vals, bode_vals, "sprial_inductor_dpa", "Spiral Inductor - DPA estimation", legend_list)
%% DPA
[lambda, x, y] = dpa(E, A, B, B, 0.8+100i, 1e-5);
Ri = (B'*x)*(y'*(E\B));
H_hat_dpa = @(s) Ri./(s-lambda);
bode_vals_dpa = H_hat_dpa(s_vals);
legend_list = [legend_list, "$\hat{H}_{\texttt{DPA}}(s)$"];
plot_bode(omega_vals, bode_vals_dpa, "sprial_inductor_dpa", "Spiral Inductor - DPA estimation", legend_list, false)

%% SADPA

opt = struct;
opt.use_lu = 0;
opt.dpa_bordered = 0;
opt.strategy = 'LM';
opt.nwanted = 20;
%opt.kmin = opt.nwanted; opt.kmax=opt.nwanted+15;
[poles, Ri, xs, ys,nr_solves]=sadpa(A, E, B, B, 0, 0.8+1i, opt);
used = [1:length(poles)];
H_hat_sadpa = @(s) sum(Ri(used)./(s-poles(used)));
bode_vals_sadpa = H_hat_sadpa(s_vals);
legend_list = [legend_list, "$\hat{H}_{\texttt{SADPA}}(s)$"];
plot_bode(omega_vals, bode_vals_sadpa, "sprial_inductor_dpa_sadpa", "Spiral Inductor - DPA estimation", legend_list, false)
fprintf("\\texttt{SADPA} usin %i poles, which it estimated using %i LU solves\n", length(poles), nr_solves)
opts = struct2table(opt);
table2latex(opts, './tex/tables/sadpa_opt.tex')


%% New start 
N=50;
% 1 
omega_exp_range = [0, 10];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
[A, B, E] = read_system(1);
H_og = @(s) B'*((E*s-A)\B);
bode_vals = bode_from_system(A, E, B, B, s_vals);
legend_list = ["$H(s)$"];
plot_bode(omega_vals, bode_vals, "sprial_inductor_irka", "Spiral Inductor - IRKA estimation", legend_list)

%% IRKA
bode_vals_irka = zeros(N, 0);
ns = [1, 2, 5, 7];
for i = 1:length(ns)
    n = ns(i);
    [Ahat, Ehat, bhat, chat, V] = irka(A, E, B, B, rand(n,1)+1i, 1e-6);
    bode_vals_irka(:, i) = bode_from_system(Ahat, Ehat, bhat, chat, s_vals);
    legend_list = [legend_list, sprintf("$\\hat{H}_{\\texttt{IRKA}}(s)$ ($n$=%i)", n)];
    plot_bode(omega_vals, bode_vals_irka(:, i), "sprial_inductor_irka", "Spiral Inductor - IRKA estimation", legend_list, false)
end

%% GRKA
[Ahat, Ehat, bhat, chat, i_val] = grka(A, E, B, B, -1+1i, 1e4, 1e8, 20, 1e-6);
bode_vals_grka = bode_from_system(Ahat, Ehat, bhat, chat, s_vals);
legend_list = [legend_list, sprintf("$\\hat{H}_{\\texttt{GRKA}}(s)$ ($i$=%i)", i_val)];
plot_bode(omega_vals, bode_vals_grka, "sprial_inductor_irka", "Spiral Inductor - IRKA estimation", legend_list, false)
