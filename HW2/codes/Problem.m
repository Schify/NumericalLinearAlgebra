classdef Problem
    %PROBLEM definition and data stored in one
    properties
        method
        name
        func
        parname="\lambda"

        A
        b
        n %size of x
        m %size of b
        p %size(s,1)

        U
        s
        sm
        V
        L = []
        X = []

        Xs
        pars %lambda or k
        rho = [] %residual norm
        eta = [] %solution norm
        F = [] %each column is one set of filter factors corresponding to lambda_i
        gcv_vals = []
        ind_best_L = []
        ind_best_gcv = []
    end
    
    methods
        function obj = Problem(A, b, method,name,L, pars, U, s, V, X)
            %constructor
            obj.method = method; obj.name = name;
            obj.A = A;
            obj.b = b;
            obj.n = size(A, 2);
            obj.m = size(A, 1);
      
            if nargin<7
                if isempty(L)
                    [U, s, V] = csvd(A);
                else
                    [U, s, X, V] = cgsvd(A, L);
                end
                
            end
            if isempty(L)
               obj.U = U; obj.s = s; obj.V = V;
            else
               obj.U = U; obj.sm = s; obj.V = V; obj.X = X;
            end
            
            obj.p = size(s,1);

            obj.L = L;
            obj.pars = pars;
           
            if strcmp(obj.method, 'Tikh')
                obj.func = @tikhonov;
            end
            if strcmp(obj.method, 'tsvd')
                obj.parname = "k";%discreete parameter
                if isempty(L)
                   obj.func = @tsvd;
                else
                   obj.func = @tgsvd;
                end
            end
            if strcmp(obj.method, 'dsvd')
                obj.func = @dsvd;
            end
            if strcmp(obj.method, 'cgls')
                if ~isempty(L)
                    error("There cannot be an L with GLS")
                end
                obj.parname = "k";%discreete parameter
                obj.pars = 1:max(obj.pars); 
                obj.func = @cgls;
            end
        end

        function obj=gen_data(obj)
            npars = length(obj.pars);
            obj.eta = zeros(npars,1);
            obj.rho = zeros(npars,1);
            obj.Xs = zeros(obj.n,npars);
            %% Solutions
            if strcmp(obj.method, 'Tikh') || strcmp(obj.method, 'tikh') || strcmp(obj.method, 'tsvd') || strcmp(obj.method, 'dsvd')
                if isempty(obj.L)%difference between general and nongeneral
                    [x_lambda, rho, eta] = obj.func(obj.U, obj.s, obj.V, obj.b, obj.pars); %#ok<PROP> 
                else
                    [x_lambda, rho, eta] = obj.func(obj.U, obj.sm, obj.X, obj.b, obj.pars); %#ok<PROP> 
                end
            end
            if strcmp(obj.method, 'cgls')
                [x_lambda, rho, eta, F] = cgls(obj.A, obj.b, max(obj.pars), 0, obj.s); %#ok<PROP> 
            end
            %% Filter factors
            if strcmp(obj.method, 'Tikh') || strcmp(obj.method, 'tsvd') || strcmp(obj.method, 'dsvd')
                if isempty(obj.L)
                    F = fil_fac(obj.s,obj.pars, obj.method); %#ok<PROP> 
                else
                    F = fil_fac(obj.sm, obj.pars, obj.method); %#ok<PROP> 
                end
            end

            obj.Xs = x_lambda;
            obj.rho = rho; %#ok<PROP> 
            obj.eta = eta; %#ok<PROP> 
            obj.F = F; %#ok<PROP> 
            % L curve
            [obj.ind_best_L,~] = find_largest_curvature(obj.rho,obj.eta, obj.pars);
            % GCV
            if strcmp(obj.method, 'cgls')
                obj.gcv_vals = obj.rho.^2./(obj.m-obj.n+obj.p-sum(obj.F)').^2;
            else
                if isempty(obj.L)%difference between general and nongeneral
                    [~,G,~] = gcv_mod(obj.U,obj.s,obj.b, obj.pars, obj.method);
                else
                    [~,G,~] = gcv_mod(obj.U,obj.sm,obj.b, obj.pars, obj.method);
                end
                close(gcf);
                obj.gcv_vals = G;
            end
            [~, obj.ind_best_gcv] = min(obj.gcv_vals);

        end
       
        function plot_l_curve(obj, fig, ax)
            if nargin<3
              fig = figure;
              ax = gca;
            end
            
            fig;

            grid on 
            set(ax, "YScale", "log")
            set(ax, "XScale", "log")
            mark_style = "--";
            if strcmp(obj.parname, "k")
                mark_style = "*";
            end
            curve=plot(ax, obj.rho, obj.eta, mark_style);
            best_ind = obj.ind_best_L; best_lambda = obj.pars(best_ind);
            best_text = sprintf("$%s = %g$",obj.parname, best_lambda);
            best_point = scatter(ax, obj.rho(best_ind), obj.eta(best_ind),100, "X");

            xlabel("$\left\| A \textbf{x} - \textbf{b}\right\|$","Interpreter","latex")
            ylabel("$\left\| \textbf{x} \right\|$","Interpreter","latex")
            
            
            if nargin<3
                title(sprintf("%s (best %s)", ...
                    obj.name, best_text), ...
                    "Interpreter","latex");
            else
                title("L curves");
                set(best_point, "DisplayName", best_text);
                set(curve, "DisplayName", obj.name);
                legend("Interpreter","latex")
            end
            



        end


        function plot_gcv_curve(obj, fig, ax)
            if nargin<3
              fig = figure;
              ax = gca;
            end
            
            fig;

            grid on 
            set(ax, "XScale", "log")
            set(ax, "YScale", "log")
            mark_style = "--";
            if strcmp(obj.parname, "k")
                mark_style = "*";
            end
            curve=plot(ax, obj.pars, obj.gcv_vals, mark_style);
            best_ind = obj.ind_best_gcv; best_lambda = obj.pars(best_ind);
            best_text = sprintf("$%s = %g$",obj.parname, best_lambda);
            best_point = scatter(ax, obj.pars(best_ind), obj.gcv_vals(best_ind),100, "X");

            xlabel(sprintf("$%s$", obj.parname),"Interpreter","latex")
            ylabel("$G$","Interpreter","latex")
            
            
            if nargin<3
                title(sprintf("%s (best %s)", ...
                    obj.name, best_text), ...
                    "Interpreter","latex");
            else
                title("GCV");
                set(best_point, "DisplayName", best_text);
                set(curve, "DisplayName", obj.name);
                legend("Interpreter","latex")
            end
            



        end

        function plot_best_sol(obj,ts,type_of_sol, fig, ax)
            
            fig;

            grid on 
%             set(ax, "XScale", "log")
%             set(ax, "YScale", "log")
            
            if strcmp(type_of_sol, 'GCV') || type_of_sol == 2
                best_ind = obj.ind_best_gcv;
                best_type = "GCV";
            elseif strcmp(type_of_sol, 'L') || type_of_sol == 1
                best_ind = obj.ind_best_L;
                best_type="L curve";
            end
            curve=plot(ax, ts, obj.Xs(:, best_ind));

            best_lambda = obj.pars(best_ind);
            best_text = sprintf("%s (%s)",obj.name, best_type);
            
            xlabel("$t$","Interpreter","latex")
            ylabel("$x(t)$","Interpreter","latex")
            
            
            title("Solution");
            set(curve, "DisplayName", best_text);
            legend("Interpreter","latex")
            



        end

    end
end

