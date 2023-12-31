function [Q_k,T_k,r,err_ind,w_k_inf_quasi] = Lanczos3(A,kmax,r,nrm_A)
 
if kmax == -1
    kmax=size(A,1);
end

% INIT
n=size(A,1);
eta = ((eps)^(3/4))/sqrt(kmax)/100;   % intermed. orth level
delta = sqrt(eps/kmax)/10;         % threshold \delta for semi-orthogonality
reorth_prev = 2;                % RECOMMENDED: reorth. against prev 2 vectors
                                % you can turn this off, but results can be
                                % quite bad
indicies = [];
err_ind = 0;                    % output error flag

%roundoff err. estim.
eps1 = sqrt(n)*eps/2;           % eps_mach = eps/2 !!! (in eta and delta 2*eps_mach is used)
gamma = 1/sqrt(2);		 % factor to check sufficient reduction
% Alloc. 
alpha = zeros(kmax+1,1);  beta = zeros(kmax+1,1);
Q_k = zeros(n,kmax);
q = zeros(n,1); beta(1)=norm(r);
w_k_inf_quasi = zeros(kmax,1);

full_reorth = 0;
w_old = [];
w = [1];
for j=1:kmax
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%    PART ONE: REGULAR LANCZOS     %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  q_old = q;
  if beta(j)==0
    q = r;
  else
    q = r / beta(j);
  end
  Q_k(:,j) = q;
  u = A*q;
  r = u - beta(j)*q_old;
  alpha(j) = q'*r;
  r = r - alpha(j)*q;
  beta(j+1) = norm(r);
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  PART THREE: COMPUTE ORTH. LOSS  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % In this part you should calculate loss of orth.
  
  % Lanczos1.m                  : compute w_{j,\infty}
  % Lanczos2.m & Lanczos3.m     : compute w vec. through recurrence (see 'update_w' below)

  if (j<kmax-1) && j>1
      %theoretically we only care about the last row of the matrix at any
      %given moment
      [w,w_old] = update_w(w, w_old , alpha, beta, n);
      w_k_inf_quasi(j+1)=max(abs(w(1:(end-1))));
      if w_k_inf_quasi(j+1)>delta
          fprintf("Reorthogonlaizing needed: %i\n", j)
          new_indicies = find_L(abs(Q_k(:, 1:j-1)'*r./norm(r)),abs(w), delta, eta);
          if isempty(indicies)
              indicies = new_indicies;
          else
              indicies = reshape(union(indicies, new_indicies), 1,[]);
          end
%           display(indicies)
%           figure
%           plot_Qq(Q_k(:,1:j),j,norm(r),r,eta,delta, abs(w), new_indicies)
          full_reorth = reorth_prev;
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%        PART FOUR: REORTH.        %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % In this part you should reorthogonalize whenever loss. of orth. is
  % detected
  
  % Lanczos1.m & Lanczos2.m     : full reorth (see int full_reorth)
  if full_reorth>0 && j>2 && j<kmax-1
    %reorthogonlaize q_j
    %Q_k, r,norm_r, delta, eta
    [new_r, new_beta] = reorth(Q_k(:, 1:j), r, norm(r), indicies);
    r = new_r;
    beta(j+1) = new_beta;
    
    full_reorth = full_reorth-1;
    
    %% Reducing the w-inf estimate, since supposedly the new vector is now orthogonal
    z3 = 3/2*randn(size(w,1)+2,1);
    w(indicies) = z3(indicies).*eps;
    if full_reorth ==0 
        indicies = [];
    end
  end


  
  
  
end

T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
end

function [w,w_old] = update_w(w,w_old,alpha,beta, n)
    % implement w recurrence here
    %adding initial zeros for k=0
    j = size(w,1);
    w = [0;w];
    w_old = [0;w_old];
    beta = [0; beta];
    alpha = [0; alpha];
    w_new = zeros(j+2, 1);%alloc
    %The diagonal and off-diagonal element in W:
    w_new(end) = 1;%diag
    w_new(end-1) = sqrt(n)*eps*beta(2)/beta(j+1)*(0.6*randn());%\psi_{j} on the offdiagonal
    for k = j:-1:2
        %because of the extension of the vectors the loop runs a bit
        %different
        % k here corresponds to k+1 in the report
        theta_jk = eps*(beta(j+1)+beta(k))*0.3*randn();%\theta_{j,k}
        w_new(k) = (beta(k)*w(k+1)+ ...
                    (alpha(k)-alpha(j+1))*w(k)+ ...
                    +beta(k-1)*w(k-1)...
                    -beta(j)*w_old(k)...
                    + theta_jk)...
                        /beta(j+1);
    end
    
    %dropping off the virtual zeros
    w_old = w(2:end);
    w     = w_new(2:end);
end

function [rnew, nrmnew, indicies] = reorth(Q_k, r,norm_r, indicies)
    Qr = Q_k(:,indicies)'*r;
%     indicies = find_L(Qq,w, delta, eta);

    if isempty(indicies)
        alpha_loc = 0;
    else
        r = r - Q_k(:,indicies)*Qr;
        alpha_loc = Q_k(:,end)'*r;
    end
    rnew = r - alpha_loc*Q_k(:,end);
    nrmnew = norm(rnew);
end

function indicies = find_L(x, w, delta, eta)
    if eta >= delta
        error('eta should be less than delta');
    end
    x = abs(x);
    w = abs(w(1:length(x)));
    above_eta   = w > eta;
    %finding contiguous intervals where eta<x
    interval_start = find(and(not([false; above_eta(1:end-1)]),...
                                above_eta)); %where could an interval start (before no, now yes)
    interval_end = find(and(above_eta,...
                         not([above_eta(2:end);false]))); %where could an interval end (now yes, next no)
    interval_needed = false(size(interval_start));
    for ind = find(x>delta)'
        interval_needed = or(interval_needed, and(ind>=interval_start, ind<=interval_end));
    end
    indicies = [];
    for int_ind = find(interval_needed)'
        indicies = [indicies, (interval_start(int_ind):interval_end(int_ind))];
    end
end

