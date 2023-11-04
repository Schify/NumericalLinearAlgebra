function plot_W(W_k, w_k_inf, x0, y0, width, height)
    subplot(1,2,1)
    ax=imagesc(log10(abs(W_k)));
    
    colormap("jet")
    title("$W_k$ : $\log_{10} {\left| {\left[W_k - I\right]}_{(i,j)}\right|}$", Interpreter="latex")
    xlabel("$i$", Interpreter="latex")
    ylabel("number of eigenvalues in interval", Interpreter="latex")
    colorbar
    set(gcf,'position',[x0,y0,width,height])
    axis equal
    
    subplot(1,2,2)
    title("$W_{k, \infty}$", Interpreter="latex")
    semilogy(w_k_inf)
    xlabel("$k$", Interpreter="latex")
    ylabel("$w_{k, \infty}$", Interpreter="latex")
    grid on

end