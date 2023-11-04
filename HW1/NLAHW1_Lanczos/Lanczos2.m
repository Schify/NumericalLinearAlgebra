function [Q_k,T_k,r,err_ind] = Lanczos2(A,kmax,r,nrm_A)
 

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
full_reorth=1;			%toggle when needed
W_k = zeros(kmax,kmax);
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

  if (j<kmax-1)
      w_next = update_w(w, w_old, alpha, beta)
      if w_k_inf(j+1)>delta
          fprintf("Reorthogonlaizing needed: %i\n", j)
          full_reorth = 2;
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%        PART FOUR: REORTH.        %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % In this part you should reorthogonalize whenever loss. of orth. is
  % detected
  
  % Lanczos1.m & Lanczos2.m     : full reorth (see boolean full_reorth)
  % Lanczos3.m                  : select L and orth w.r.t. vectors in L
  %                               DON'T FORGET: update w afterwards!
  % DON'T FORGET: in next run of main loop, reorth against L (or full) again!!!
  % Write your reorth in a function of the form 
  % [r,nrmnew] = reorth(Q_k,r,nrm,L), with nrm the norm of r on input, nrmnew norm on output;
  
  
  
  
  
end

T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
end

function [w,w_old] = update_w(w,w_old,alpha,beta)
    % implement w recurrence here
end


