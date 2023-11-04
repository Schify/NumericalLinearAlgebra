function [Q_k,T_k,r,err_ind,w_k_inf_quasi] = Lanczos2(A,kmax,r,nrm_A)
 

% INIT
n=size(A,1);
eta = ((eps)^(3/4))/sqrt(kmax);   % intermed. orth level
delta = sqrt(eps/kmax);         % threshold \delta for semi-orthogonality
reorth_prev = 1;                % RECOMMENDED: reorth. against prev 2 vectors
                                % you can turn this off, but results can be
                                % quite bad
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
      if j == 2
          w_old = [];
          w = [1];
      end
      [w,w_old] = update_w(w, w_old , alpha, beta, n);
      w_k_inf_quasi(j)=max(abs(w(1:(end-1))));
      if w_k_inf_quasi(j)>delta
          fprintf("Reorthogonlaizing needed: %i\n", j)
          where_update = abs(w(1:(end-1)))>delta;
          z3 = 3/2*randn(size(w));
          w(where_update) = z3(where_update)*eps;

          full_reorth = 2;
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%        PART FOUR: REORTH.        %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % In this part you should reorthogonalize whenever loss. of orth. is
  % detected
  
  % Lanczos1.m & Lanczos2.m     : full reorth (see int full_reorth)
  if full_reorth>0 && j>2
    %reorthogonlaize q_j
    fprintf("Reorth")
    display(r(1:5))
    display(norm(Q_k(:,1:j-1)'*r)/j)

    [new_r, new_beta] = reorth(Q_k(:, 1:j), r, norm(r));
    r = new_r;
    beta(j+1) = new_beta;

    full_reorth = full_reorth-1;

    display(r(1:5))
    display(norm(Q_k(:,1:j-1)'*r)/j)
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
    w_new = zeros(j+2, 1);%alloc
    
    %The diagonal and off-diagonal element in W:
    w_new(end) = 1;%diag
    w_new(end-1) = sqrt(n)*eps*beta(1)/beta(j)*(0.6*randn());%\psi_{j} on the offdiagonal
    for k = j:-1:2
        %because of the extension of the vectors the loop runs a bit
        %different
        % k here corresponds to k+1 in the report
        theta_jk = eps*(beta(j)+beta(k-1))*0.3*randn();%\theta_{j,k}
        w_new(k) = (beta(k-1)*w(k+1)+(alpha(k-1)-alpha(j))*w(k)+ ...
                    -beta(j-1)*w_old(k) + theta_jk)...
                        /beta(j);
    end
    
    %dropping off the virtual zeros
    w_old = w(2:end);
    w     = w_new(2:end);

end

function [rnew, nrmnew] = reorth(Q_k, r, nrm)
    r = r - Q_k(:,1:(end-1))*(Q_k(:,1:(end-1))'*r);
    %r = r - Q_k(:, 1:j-1)*(Q_k(:, 1:j-1)'*r);
    alpha_loc = Q_k(:,end)'*r;
    rnew = r - alpha_loc*Q_k(:,end);
    nrmnew = norm(rnew);
end

