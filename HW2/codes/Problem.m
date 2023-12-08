classdef Problem
    %PROBLEM definition and data stored in one
    properties
        method
        name
        func

        A
        b
        n %size of x
        m %size of b

        U
        s
        sm
        V
        L = []

        X = []
        pars %lambda or k
        rho = [] %residual norm
        eta = [] %solution norm
        F = [] %each column is one set of filter factors corresponding to lambda_i
    end
    
    methods
        function obj = Problem(A, b, method,name,L, pars, U, s, V, X)
            %constructor
            obj.method = method; obj.name = name;
            obj.A = A;
            obj.b = b;
            obj.n = size(A, 2);
            obj.m = size(A, 1);
      
            if nargin<6
                if isempty(L)
                    [U, s, V] = csvd(A);
                else
                    [U, sm, X, V] = csvd(A);
                end
                
            end
            if isempty(L)
               obj.U = U; obj.s = s; obj.V = V;
            else
               obj.U = U; obj.sm = sm; obj.V = V; obj.X = X;
            end
            
            obj.pars = pars;
           
            if strcmp(obj.method, 'Tikh')
                obj.func = @tikhonov;
            end
        end

        function obj=gen_data(obj)
            npars = length(obj.pars);
            obj.eta = zeros(npars);
            obj.rho = zeros(npars);
            obj.X = zeros(obj.n,npars);

            if strcmp(obj.method, 'Tikh')
                x0 = zeros(obj.n,1);
                if isempty(obj.L)
                    [x_lambda, rho, eta] = tikhonov(obj.U, obj.s, obj.V, obj.b, obj.pars, x0); %#ok<PROP> 
                else
                    [x_lambda, rho, eta] = tikhonov(obj.U, obj.sm, obj.X, obj.b, obj.pars, x0); %#ok<PROP> 
                end
            end

            obj.X = x_lambda;
            obj.rho = rho; %#ok<PROP> 
            obj.eta = eta; %#ok<PROP> 
        end
       
        function plot_l_curve(obj, fig, ax)
            put_legend = false;
            if nargin<3
              fig = figure;
              ax = gca;
              put_legend = true;
            end
            [best_ind, best_lambda] = find_largest_curvature(obj.eta, obj.rho, obj.pars);
            fig;
            hold on 
            
            grid on 
            set(ax, "XScale", "log")
            set(ax, "YScale", "log")
            curve=plot(ax, obj.eta, obj.rho);
            xlabel("$\left\| A \textbf{x} - \textbf{b}\right\|$","Interpreter","latex")
            ylabel("$\left\| \textbf{x} \right\|$","Interpreter","latex")
            scatter(ax, obj.eta(best_ind), obj.rho(best_ind))
            if put_legend
                set(curve, "DisplayName", obj.name)
                legend("Interpreter","latex")
            else
                title(obj.name,"Interpreter","latex")
            end
            

            hold off


        end

    end
end

