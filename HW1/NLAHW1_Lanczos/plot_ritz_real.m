function plot_ritz_real(T, real_eig_values, leave_out_last)
%% Plot the real part of the ritz and eigenvalues
    if nargin < 3
        leave_out_last=40;
    end
    
    
    set(gcf,'position',[200,100,850,650])
    clf
    hold on
    yyaxis left
    
    num_bins = floor(sqrt(size(real_eig_values,1))*2);
%     bin_edges = 2.^linspace(log2(min(real_eig_values)*0.9), ...
%                          log2(max(real_eig_values)*1.1), ...
%                             num_bins);
    bin_edges = linspace((min(real_eig_values)*0.9), ...
                         (max(real_eig_values)*1.1), ...
                            num_bins);
    [N,bin_edges] = histcounts(real(real_eig_values),bin_edges);
    bar((bin_edges(2:end)+bin_edges(1:end-1))/2.0,N,0.9);
   
    xlabel("$\mathcal{R}\left(\lambda_i\right)$ or $\mathcal{R}\left(\mu_i\right)$", Interpreter="latex")
    ylabel("number of eigenvalues ($\lambda_i$)", Interpreter="latex")
    ylim([0.1 inf])
    set(gca, "YScale", "log")
    
    yyaxis right
    x_vals = [];
    y_vals = [];
    for j = 1:(size(T,1)-leave_out_last)
        fprintf("%i \n", j)
        Diag_Tk = eig(T(1:j, 1:j));
        y_vals = [y_vals; j*ones(j,1)];
        x_vals = [x_vals; real(Diag_Tk)];
    end %for
    scatter(x_vals,y_vals ,"rx")
    ylabel("j", Interpreter="latex") 
    
    writematrix([x_vals, y_vals], "../plots/Lanczos0_ritz.csv")
    writematrix([((bin_edges(2:end)+bin_edges(1:end-1))/2.0)',N'], ...
        "../plots/Lanczos0_lambda.csv")
    %set(gca, "XScale", "log")
end

