function [outputArg1,outputArg2] = plot_Qq(Q_k,j,norm_r,r,eta, delta, w, indicies) 
    plot(abs(Q_k(:,1:j)'*r./norm_r))
    hold on 
    plot(abs(w))
    plot(indicies, w(indicies), "x", "Color","magenta")
    legend(["$Q_{j}^*\mathbf{q}_{j+1}$", "$w$", "selected"], Interpreter="latex", AutoUpdate='off', Location="southeast")
    
    hold off
    title(['$\left|Q_j^* \mathbf{q}_{j+1}\right|$ with $j=$', num2str(j)], Interpreter="latex")
    xlabel("Indicies")
    set(gca, "YScale", "log")
    grid on
%     xline(reorth_ind)
%     xline(reorth_ind, "r--")

    yline(eta, "g--", "$\eta$", "Alpha",0.7,"Interpreter","latex")
    yline(delta, "o--", "$\delta$", "Alpha",0.7,"Interpreter","latex")
end

