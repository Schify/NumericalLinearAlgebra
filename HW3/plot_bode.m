function [outputArg1,outputArg2] = plot_bode(omega, bode_vals, name, full_name, legend_names,newfig)
%PLOT_BODE create and save bode plot for data
if nargin < 6 || newfig
    figure
    plot_new = true;
    linetype = "-";
else 
    hold on
    plot_new = false;
    linetype = "--";
end 

subplot(2,1,1)
if plot_new
    loglog(omega, abs(bode_vals),linetype)
else
    hold on
    plot(omega, abs(bode_vals),linetype)
end
subtitle("Amplitude")
ylabel("$\left|H(2\pi\omega i)\right|$", Interpreter="latex")
grid on

title(sprintf("Bode Plot - %s", full_name))

subplot(2,1,2)

semilogx(omega, angle(bode_vals),linetype)
subtitle("Phase")
ylabel("$\arg\left(H(2\pi\omega i)\right)$", Interpreter="latex")
xlabel("$\omega$", Interpreter="latex")
grid on

if min(size(bode_vals))>1 || (nargin >4 && ~isempty(legend_names))
    if nargin <= 4 || isempty(legend_names)
        legend_names = "out ch"+[1:min(size(bode_vals))];
    end
    for i=2:2
        subplot(2,1,i)
        legend(legend_names,"Location","northwest", "Interpreter","latex")
    end
end

matlab2tikz(sprintf('./tex/figs/bode_%s.tex', name), 'extraTikzpictureOptions', ['thick, trim axis left , trim axis right'],...
    'width','0.95\textwidth')

hold off

end

