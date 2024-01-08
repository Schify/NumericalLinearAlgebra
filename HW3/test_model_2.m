clear 
close all 
addpath("../matlab2tikz/src/")
addpath("./methods/")
addpath("../")

N=50;
output_ch = 2;
%% 2
omega_exp_range = [0, 4];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
beta = 1e-6;
[B, C, K, M, C_names ] = read_system(2);
C = C(output_ch,:)'; C_names = C_names(output_ch);%restrict to one output channel
D = beta*K;
foo = @(s) M*s^2 + D*s+K;
bode_vals = bode_from_function(foo,B,C, s_vals)';
legend_list = ["$H(s)$"];
plot_bode(omega_vals, bode_vals, "butterfly_gyro_sadpa", "Butterfly Gyroscope", C_names)

% %% system conversion
% M_size = size(M);
% E = [eye(M_size), zeros(M_size); zeros(M_size), M];
% A = [zeros(M_size), eye(M_size); -K, -K];
% B_ext = [zeros(M_size(1), 1); B];
% C_ext = [zeros(size(C)), C]';

% %% QDPA
% [lambda, x, y] = qdpa(M, D, K, B, C, -0.1+1i, 1e-4);
% Ri = (C'*x)*(y'*(B));
% H_hat_dpa = @(s) Ri./(s-lambda);
% bode_vals_dpa = H_hat_dpa(s_vals);
% legend_list = [legend_list, "$\hat{H}_{\texttt{QDPA}}(s)$"];
% plot_bode(omega_vals, bode_vals_dpa, "gyro_qdpa", "Butterfly Gyroscope - QDPA estimation", legend_list, false)

% %% QDPA
% 
% opt = struct;
% opt.use_lu = 0;
% opt.dpa_bordered = 0;
% opt.strategy = 'LM';
% opt.nwanted = 5;
% %opt.kmin = opt.nwanted; opt.kmax=opt.nwanted+15;
% [poles, Ri, xs, ys,nr_solves]=sadpa(A, E,  B_ext, C_ext, 0, 0.8+1i, opt);
% used = [1:length(poles)];
% H_hat_sadpa = @(s) sum(Ri(used)./(s-poles(used)));
% bode_vals_sadpa = H_hat_sadpa(s_vals);
% legend_list = [legend_list, "$\hat{H}_{\texttt{SADPA}}(s)$"];
% plot_bode(omega_vals, bode_vals_sadpa, "butterfly_gyro_sadpa", "Butterfly Gyroscope - SADPA estimation", legend_list, false)
% fprintf("\\texttt{SADPA} usin %i poles, which it estimated using %i LU solves\n", length(poles), nr_solves)
% opts = struct2table(opt);
% table2latex(opts, './tex/tables/sadpa_butterfly_gyro_opt.tex')

%% IRKA
Ahat, Ehat, bhat, chat, V] = irka(M, D, K, B, C, -0.1+1i, 1e-4);
Ri = (C'*x)*(y'*(B));
H_hat_dpa = @(s) Ri./(s-lambda);
bode_vals_dpa = H_hat_dpa(s_vals);
legend_list = [legend_list, "$\hat{H}_{\texttt{QDPA}}(s)$"];
plot_bode(omega_vals, bode_vals_dpa, "gyro_qdpa", "Butterfly Gyroscope - QDPA estimation", legend_list, false)

