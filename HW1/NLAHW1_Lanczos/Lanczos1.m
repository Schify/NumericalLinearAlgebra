function [Q_k,T_k,r,err_ind, w_k_inf] = Lanczos0(A,kmax,r,nrm_A)
%Theoretical Lanczos algortihm
    calc_w_k_inf = false;
    if nargout > 4
        calc_w_k_inf = true;
    end


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
w_k_inf = zeros(kmax, 1);

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
  if (j<kmax-1)
      if beta(j+1)==0
        q_new = r;
      else
        q_new = r / beta(j+1);
      end
      Q_k(:, j+1) = q_new;
      w_k_inf(j+1) = max(abs(Q_k(:, 1:j+1)'*Q_k(:, 1:j+1)-eye(j+1,j+1)), [], "all");
      if w_k_inf(j+1)>delta
          fprintf("Reorthogonlaizing needed: %i\n", j)
          full_reorth = 2;
      end
  end
  % Lanczos2.m & Lanczos3.m     : compute w vec. through recurrence (see 'update_w' below)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%        PART FOUR: REORTH.        %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % In this part you should reorthogonalize whenever loss. of orth. is
  % detected
  
  % Lanczos1.m & Lanczos2.m     : full reorth (see boolean full_reorth)
  if full_reorth>0 && j>2
    %reorthogonlaize q_j
    r = r - Q_k(:, 1:j-1)*(Q_k(:, 1:j-1)'*r);
    %r = r - Q_k(:, 1:j-1)*(Q_k(:, 1:j-1)'*r);
    alpha_loc = q'*r;
    r = r - alpha_loc*q;
    beta(j+1) = norm(r);
    full_reorth = full_reorth-1;
  end
  

 end

T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
end


