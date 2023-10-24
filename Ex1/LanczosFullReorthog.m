% Matlab code LanczosFullReorthog.m
% For "Applied Numerical Linear Algebra",  Chapter 7, section 3
% Written by James Demmel, Jun  6, 1997
%
% Perform Lanczos with complete reorthogonalization
%
% Inputs:
%
%   A = input matrix (may be sparse)
%   e = true eigenvalues, sorted from largest to smallest
%   q = starting vector for Lanczos
%   m = number of Lanczos steps to take
%   mm1,mm2 = indices of extreme eigenvalues for which to plot 
%             error and error bounds
%   mine, maxe = range in which to plot eigenvalues
%   steplabel = array of values at which to draw vertical lines in plots
%   recomp = 1 to perform algorithm from scratch
%          = 0 just to redisplay data, with possibly different
%               mine, maxe, mm1, mm2, steplabel
%   whichbound = 2 to plot all errors and error bound
%              = 1 to plot local error and error bounds
%              = 0 to plot error bounds only
%
% Outputs:
%
%   plot of Ritz values for all m Lanczos steps, in range [mine,maxe]
%
%   plot of global error (distance from final eigenvalue) (if whichbound=2)
%           local error (distance from closest eigenvalue) (if whichbound=1)
%           error bound (provided by Lanczos algorithm)
%     for mm1 to mm2 smallest and mm1 to mm2 largest eigenvalues
%
%   plot of components of each eigenvector in current Lanczos vector
%           (only valid if A is diagonal with eigenvalues sorted
%            from largest to smallest!)
%     for mm1 to mm2 smallest and mm1 to mm2 largest eigenvalues
%
linewidth = 2;
linewidthd = 6;
markersize = 8;
if (recomp==1)  % {
   n=max(size(A));
   Eval = zeros(m,m);
   ErrorBounds = zeros(m,m);
   qq = q/norm(q);
   Q = qq;
   clear alpha
   clear beta

   disp('start performing Lanczos with full reorthogonalization')
   disp('print step number')
   for i=1:m
      if (rem(i,10)==0)
         i,
      end
      z = A*qq;
      alpha(i) = qq'*z;
%     reorthogonalize twice to guarantee orthogonality of Lanczos vectors
      z = z - Q*(Q'*z);
      z = z - Q*(Q'*z);
      beta(i) = norm(z);
      qq = z/beta(i);
      Q = [Q,qq];
   end

   disp('start extracting Ritz values, print step number')
   T = diag(alpha(1:m)) + diag(beta(1:m-1),1) + diag(beta(1:m-1),-1);
   for i=1:m
     if (rem(i,10)==0)
        i,
     end
     [V,D] = eig(T(1:i,1:i));
%    sort eigenvalues into decreasing order
     [Ds,Is] = sort(-diag(D));
     Eval(1:i,i) = -Ds;
     Vs = V(:,Is);
     ErrorBounds(1:i,i) = abs((beta(i)*Vs(i,:))');
   end
   disp('largest deviation of singular value of Lanczos vector matrix Q from 1 ='),
      max(abs(svd(Q)-1))
   disp('done with extracting eigenvalues, bounds')
   disp('pause (hit return to continue)'), pause
end % }
%
% Produce graph of converging Ritz values 
%
figure(1)
hold off, clf
for j=1:ceil(m/2),
    rr = 2*j-1;
    if (rem(j,4)==1),
       str='k+';
    elseif (rem(j,4)==2),
       str='r+';
    elseif (rem(j,4)==3),
       str='g+';
    elseif (rem(j,4)==0),
       str='b+';
    end
    deval = diag(Eval,j-1);deval=deval(j:m-j+1);
    hndl=plot(rr:m,Eval(j,rr:m),str); 
    set(hndl,'MarkerSize',markersize);
    hold on
    hndl=plot(rr:m,deval,str);
    set(hndl,'MarkerSize',markersize);
end
hndl=plot((m+1)*ones(size(e)),e,'k+');
set(hndl,'MarkerSize',markersize);
axis([0,m+2,mine,maxe])
for b=steplabel, plot([b,b],[mine,maxe],'k:'); end
title([int2str(m),' steps of Lanczos (full reorthogonalization) applied to A'])
xlabel('Lanczos step'), ylabel('Eigenvalues')
disp('done with figure 1')
disp('pause (hit return to continue)'),pause
%
% Produce graph of errors, error bounds in largest mm1 to mm2 Ritz values
%
figure(2)
hold off, clf
step = 3;
for js=mm1:step+1:mm2
  hold off, clf
  for j=js:min(js+step,mm2)
    if (rem(j,4)==1),
       str='k';
       strB='--k';
       strD=':k';
    elseif (rem(j,4)==2),
       str='r';
       strB='--r';
       strD=':r';
    elseif (rem(j,4)==3),
       str='g';
       strB='--g';
       strD=':g';
    elseif (rem(j,4)==0),
       str='b';
       strB='--b';
       strD=':b';
    end
    disp(['Eigenvalue being plotted = ',num2str(e(j)), ' :: color = ', str]),
%   plot global error
    if (whichbound > 1),
       errpl = abs(Eval(j,j:m)-e(j))/(abs(e(j)));
       errpl = max( errpl, 1e-17 );
       hndl=semilogy((j:m),errpl,str);
       set(hndl,'LineWidth',linewidth);
       set(hndl,'MarkerSize',markersize);
       hold on
    end
%   plot local error
    if (whichbound > 0),
       errpl = min(abs(ones(size(e))*Eval(j,j:m)-e*ones(1,m-j+1)))/(abs(e(j)));
       errpl = max( errpl, 1e-17 );
       hndl=semilogy((j:m), errpl, strD );
       set(hndl,'LineWidth',linewidthd);
       set(hndl,'MarkerSize',markersize);
       hold on
    end
%   plot error bound
    errpl = ErrorBounds(j,j:m)/(abs(e(j)));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy((j:m),errpl,strB);
    set(hndl,'LineWidth',linewidth);
    set(hndl,'MarkerSize',markersize);
    hold on
  end
  axis([0,m+1,1e-16,1])
  jsend = min(js+step,mm2);
  if (whichbound > 1),
   if (jsend ~= js),
      title(['True errors and error bounds, in eigenvalues ',int2str(js), ...
          ' to ', int2str(min(js+step,mm2))])
   else
      title(['True errors and error bounds, in eigenvalue ',int2str(js)])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
   ylabel('Global Error (solid), Local Error (dotted), Bound (dashed)')
  elseif (whichbound > 0),
   if (jsend ~= js),
      title(['True local error and error bounds, in eigenvalues ',int2str(js), ...
          ' to ', int2str(min(js+step,mm2))])
   else
      title(['True local error and error bounds, in eigenvalue ',int2str(js)])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
   ylabel('Local Error (dotted), Bound (dashed)')
  else
   if (jsend ~= js),
      title(['Error bounds, in eigenvalues ',int2str(js), ...
          ' to ', int2str(min(js+step,mm2))])
   else
      title(['Error bounds, in eigenvalue ',int2str(js)])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
  end
  for b=steplabel, plot([b,b],[1e-16,1],'k:'); end
  disp('pause (hit return to continue)'), pause
end
disp('done with figure 2')
disp('pause (hit return to continue)'), pause
%
% Produce graph of errors, error bounds in smallest mm1 to mm2 Ritz values
%
figure(3)
hold off, clf
for js=mm1:step+1:mm2
  hold off, clf
  for j=js:min(js+step,mm2)
    if (rem(j,4)==1),
       str='k';
       strB='--k';
       strD=':k';
    elseif (rem(j,4)==2),
       str='r';
       strB='--r';
       strD=':r';
    elseif (rem(j,4)==3),
       str='g';
       strB='--g';
       strD=':g';
    elseif (rem(j,4)==0),
       str='b';
       strB='--b';
       strD=':b';
    end
    deval = diag(Eval,j-1); deval=deval(1:m-j+1);
    deval = deval';
    disp(['Eigenvalue being plotted = ',num2str(e(n+1-j)), ' :: color = ', str]),
%   plot global error
    if (whichbound > 1),
       errpl = abs(deval-e(n+1-j))/max(abs(e));
       errpl = max( errpl, 1e-17 );
       hndl=semilogy((j:m), errpl,str);
       set(hndl,'LineWidth',linewidth);
       set(hndl,'MarkerSize',markersize);
       hold on
    end
%   plot local error
    if (whichbound > 0),
       errpl = min(abs(ones(size(e))*deval-e*ones(1,m-j+1)))/(abs(e(j)));
       errpl = max( errpl, 1e-17 );
       hndl=semilogy((j:m), errpl, strD );
       set(hndl,'LineWidth',linewidthd);
       set(hndl,'MarkerSize',markersize);
       hold on
    end
    devalbnd = diag(ErrorBounds,j-1); devalbnd = devalbnd(1:m-j+1);
%   plot error bound
    errpl = devalbnd/max(abs(e));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy((j:m),errpl,strB);
    set(hndl,'LineWidth',linewidth);
    set(hndl,'MarkerSize',markersize);
    hold on
  end
  for b=steplabel, plot([b,b],[1e-16,1],'k:'); end
  axis([0,m+1,1e-16,1])
  jsdif = (n+1-js)-(n+1-min(js+step,mm2));
  if (whichbound > 1),
   if (jsdif ~= 0)
      title(['True errors and error bounds, in eigenvalues ', ...
          int2str(n+1-min(js+step,mm2)),' to ', int2str(n+1-js)])
   else
      title(['True errors and error bounds, in eigenvalue ',  ...
          int2str(n+1-min(js+step,mm2))])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
   ylabel('Global Error (solid), Local Error (dotted), Bound (dashed)')
  elseif (whichbound > 0),
   if (jsdif ~= 0)
      title(['True local error and error bounds, in eigenvalues ', ...
          int2str(n+1-min(js+step,mm2)), ' to ', int2str(n+1-js)])
   else
      title(['True local error and error bounds, in eigenvalue ', ...
          int2str(n+1-min(js+step,mm2))])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
   ylabel('Local Error (dotted), Bound (dashed)')
  else
   if (jsdif ~= 0)
      title(['Error bounds, in eigenvalues ',int2str(n+1-min(js+step,mm2)), ...
          ' to ', int2str(n+1-js)])
   else
      title(['Error bounds, in eigenvalue ',int2str(n+1-min(js+step,mm2))])
   end
   xlabel('Lanczos step (full reorthogonalization)'),
  end
  disp('pause (hit return to continue)'), pause
end
disp('done with figure 3')
disp('pause (hit return to continue)'), pause
%
% Plot components of Lanczos vectors for largest mm1 to mm2 eigenvectors
%
figure(4)
hold off, clf
for js=mm1:step+1:mm2
  hold off, clf
  for j=js:min(js+step,mm2)
    if (rem(j,4)==1),
       str='k';
       strB='--k';
       strQp='+k';
       strQn='ok';
    elseif (rem(j,4)==2),
       str='r';
       strB='--r';
       strQp='+r';
       strQn='or';
    elseif (rem(j,4)==3),
       str='g';
       strB='--g';
       strQp='+g';
       strQn='og';
    elseif (rem(j,4)==0),
       str='b';
       strB='--b';
       strQp='+b';
       strQn='ob';
    end
    disp(['Eigenvalue being plotted = ',num2str(e(j)), ' :: color = ', str]),
    Qpos = find(Q(j,1:m)>0);
    Qneg = find(Q(j,1:m)<0);
    errpl = abs(Q(j,Qpos));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy(Qpos,errpl,strQp);
    set(hndl,'LineWidth',linewidth);
    set(hndl,'MarkerSize',markersize);
    hold on
    errpl = abs(Q(j,Qneg));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy(Qneg,errpl,strQn);
    set(hndl,'LineWidth',linewidth);
    set(hndl,'MarkerSize',markersize);
    axis([0,m+1,1e-16,1])
  end
  axis([0,m+1,1e-16,1])
  if (js ~= min(js+step,mm2))
     title(['Lanczos vector components for eigenvectors  ',int2str(js), ...
         ' to ', int2str(min(js+step,mm2))])
  else
     title(['Lanczos vector components for eigenvector  ',int2str(js)])
  end
  xlabel('Lanczos step (full reorthogonalization)'), ylabel('eigencomponent (+ = pos, o = neg)')
  for b=steplabel, plot([b,b],[1e-16,1],'k:'); end
  for b = [1e-12 1e-8 1e-4], plot([0,m+1],[b,b],'k:'); end
  disp('pause (hit return to continue)'), pause
end
disp('done with figure 4')
disp('pause (hit return to continue)'), pause
%
% Plot components of Lanczos vectors for smallest mm1 to mm2 eigenvectors
%
figure(5)
hold off, clf
for js=mm1:step+1:mm2
  hold off, clf
  for j=js:min(js+step,mm2)
    if (rem(j,4)==1),
       str='k';
       strB='--k';
       strQp='+k';
       strQn='ok';
    elseif (rem(j,4)==2),
       str='r';
       strB='--r';
       strQp='+r';
       strQn='or';
    elseif (rem(j,4)==3),
       str='g';
       strB='--g';
       strQp='+g';
       strQn='og';
    elseif (rem(j,4)==0),
       str='b';
       strB='--b';
       strQp='+b';
       strQn='ob';
    end
    disp(['Eigenvalue being plotted = ',num2str(e(n+1-j)), ' :: color = ', str]),
    Qpos = find(Q(n+1-j,1:m)>0);
    Qneg = find(Q(n+1-j,1:m)<0);
    errpl = abs(Q(n+1-j,Qpos));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy(Qpos,errpl,strQp);
    set(hndl,'MarkerSize',markersize);
    set(hndl,'LineWidth',linewidth);
    hold on
    errpl = abs(Q(n+1-j,Qneg));
    errpl = max( errpl, 1e-17 );
    hndl=semilogy(Qneg,errpl,strQn);
    set(hndl,'MarkerSize',markersize);
    set(hndl,'LineWidth',linewidth);
    end
  axis([0,m+1,1e-16,1])
  if (n+1-min(js+step,mm2) ~= n+1-js)
     title(['Lanczos vector components for eigenvectors ', ...
         int2str(n+1-min(js+step,mm2)), ' to ', int2str(n+1-js)])
  else
     title(['Lanczos vector components for eigenvector ',int2str(n+1-min(js+step,mm2))])
  end
  xlabel('Lanczos step (full reorthogonalization)'), ylabel('eigencomponent (+ = pos, o = neg)')
  for b=steplabel, plot([b,b],[1e-16,1],'k:'); end
  for b = [1e-12  1e-8  1e-4 ], plot([0,m+1],[b,b],'k:'); end
  disp('pause (hit return to continue)'), pause
end
disp('done with figure 5')
