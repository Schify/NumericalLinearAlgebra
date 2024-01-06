function [outputArg1,outputArg2] = plot_bode(omega, bode_vals, name, full_name, legend_names)
%PLOT_BODE create and save bode plot for data
figure

subplot(2,1,1)
loglog(omega, abs(bode_vals))
subtitle("Amplitude")
ylabel("$\left|H(2\pi\omega i)\right|$", Interpreter="latex")
grid on

title(sprintf("Bode Plot - %s", full_name))

subplot(2,1,2)

semilogx(omega, angle(bode_vals))
subtitle("Phase")
ylabel("$\arg\left(H(2\pi\omega i)\right)$", Interpreter="latex")
xlabel("$\omega$", Interpreter="latex")
grid on

if nargin > 4 || min(size(bode_vals))>1
    if nargin <= 4
        legend_names = "output_"+[1:min(size(bode_vals))];
    end
    for i=2:2
        subplot(2,1,i)
        legend(legend_names,"Location","northwest")
    end
end

matlab2tikz(sprintf('./tex/figs/bode_%s.tex', name), 'extraTikzpictureOptions', ['trim axis left , trim axis right'],...
    'width','0.95\textwidth')

end

