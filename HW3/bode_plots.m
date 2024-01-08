clear 
close all 
addpath("../matlab2tikz/src/")

N=700;
%% 1 
omega_exp_range = [0, 10];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
[A, B, E] = read_system(1);
bode_vals = bode_from_system(A, E, B, B, s_vals);
plot_bode(omega_vals, bode_vals, "sprial_inductor", "Spiral Inductor")

%% 2
omega_exp_range = [0, 4];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
beta = 1e-6;
[B, C, K, M, C_names ] = read_system(2);
foo = @(s) M*s^2 + beta*K*s+K;
bode_vals = bode_from_function(foo,B,C', s_vals)';
plot_bode(omega_vals, bode_vals, "butterfly_gyro", "Butterfly Gyroscope", C_names)

%% 3
omega_exp_range = [-2, 2];
omega_vals = logspace(omega_exp_range(1), omega_exp_range(2), N);
s_vals = 2*pi*1i*omega_vals;
beta = 1e-6;
[A, B, C, D] = read_system(3);
for input_no = 1:size(B, 2)
    bode_vals = bode_from_system(A,eye(size(A)),B(:, input_no),C(:,:)', s_vals);%for the first input
    plot_bode(omega_vals, bode_vals, sprintf("iss_i%i", input_no), sprintf("ISS A12 - Input %i",input_no))
end