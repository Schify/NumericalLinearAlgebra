function I=BF_HH_pldc_new(p,q,P,QUAD,phi,psi,kappa,i,j)
%   general routine for t-t interaction
%   does ALL WORK for the t-t interaction i.e.
%   determines #common, reordering, quad select,...
    if(i<=j)
    [common,shift,orient]=reorder_common_TOT(p',q');
    phi=mod(phi+orient,3);
    psi=mod(psi+shift,3);
    p=circshift(p,orient);
    q_temp=circshift(q,shift);
    if (common==2 && (q_temp(1)==p(2)))
    q_temp = [q_temp(2); q_temp(1); q_temp(3)];
    psi_temp = 2;
    if(psi==0)
        psi_temp = 1;
    elseif(psi==1)
        psi_temp = 0;
    end
    psi=psi_temp;
    end
    
    P0=P(p(1),:)';
    P1=P(p(2),:)';
    P2=P(p(3),:)';
    
    Q0=P(q_temp(1),:)';
    Q1=P(q_temp(2),:)';
    Q2=P(q_temp(3),:)';

    Ap=[(P1)-P0,(P2)-P1];
    Aq=[(Q1)-Q0,(Q2)-Q1];

    gp=(sqrt((det(Ap'*Ap))));
    gq=(sqrt((det(Aq'*Aq))));

    %g_vec=@(v,phi,psi)ghat(v,P0,Q0,Ap,Aq,phi,psi,kappa);
    q_SAUTER=0;
    XW=QUAD{common+1};
    c_geom=(gp*gq/(4*pi));
    for i=1:length(XW)
        g=0;
        for j = 1:length(XW(i).X)
            g=g+ghat(XW(i).X{j},P0,Q0,Ap,Aq,phi,psi,kappa);%g_vec(XW(i).X{j},phi,psi);
        end
        q_SAUTER=q_SAUTER+g*XW(i).W;
    end
    I=c_geom*q_SAUTER;
    else
        I=BF_HH_pldc_new(q,p,P,QUAD,psi,phi,kappa,j,i);
    end
end

function r = R_vec(xy,P0,Q0,Ap,Aq)
x=xy(1:2,:);
y=xy(3:4,:);
pminq=(P0+Ap*x-(Q0+Aq*y));
r=vecnorm(pminq);
end
function g=ghat(v,P0,Q0,Ap,Aq,phi,psi,kappa)
R=R_vec(v,P0,Q0,Ap,Aq);
E=exp(1i*kappa*R);
g=(E./R).*(dof_eval(v(1:2,:),phi)).*(dof_eval(v(3:4,:),psi));
end
function b=dof_eval(x,dof)
u=x(1,:);
v=x(2,:);
if(dof==0)
    b=1-u;
elseif(dof==1)
    b=u-v;
else
    b=v;
end
end