function mat=SLP_pld(H,block,assembler,kappa)
% Function that computes the Single Layer Potential Boundary operator
% w.r.t. piecewise linear DIScontinuous basis functions.
% Input:    H          -> struct representing the pld function space
%           block      -> dof clusters, contain index sets

%           assembler  ->  method defining the assembly type (e.g. 'dense','aca',...)
%                          contains quadrature as well
%                           
%           kappa      -> wave number

% STEP ONE: COMPUTE LOCAL ATTACHMENTS
mat.LR=0;
sig=block.sig;
tau=block.tau;
grid=H.grid;
P=grid.P;
tri=grid.tri;
vA=grid.vtx_attachments;
tA=grid.tri_attachments;

QUAD = assembler.QUAD;

B=assembler.B;


I=sig.indices';
N       =   numel(I);
IA=[];
for i=I
    IA=[IA,(3*(i-1)+1:3*i)];
end
J=tau.indices';
JA=[];
for j=J
    JA=[JA,(3*(j-1)+1:3*j)];
end


tA_sub=tA(I,J);
n_sig=numel(I);
n_tau=numel(J);


%N=size(P0,1)/3;
XW=QUAD{1};
X=XW(1).X{1};

P0_sig=grid.AFFINE_MAPS.P0(IA);
A_sig=grid.AFFINE_MAPS.A(IA,:);
g_sig=grid.AFFINE_MAPS.g(I);

P0_tau=grid.AFFINE_MAPS.P0(JA);
A_tau=grid.AFFINE_MAPS.A(JA,:);
g_tau=grid.AFFINE_MAPS.g(J);

XI=P0_sig+A_sig(:,1:2)*X(1:2,:);
ETA=P0_tau+A_tau(:,1:2)*X(3:4,:);

S_1m1_fast=zeros(9*n_sig,n_tau);

Bmat_TOT=[B{1,1},B{2,1},B{3,1},B{1,2},B{2,2},B{3,2},B{1,3},B{2,3},B{3,3}];

N=size(P0_sig,1)/3;
if(strcmp(assembler.method,'dense')||~block.adm)
    for j_ind=1:numel(J)
        eta=ETA(3*(j_ind-1)+1:3*j_ind,:);
        S_1m1_fast(:,j_ind)=BF_FAR(XI,eta,g_sig,g_tau(j_ind),kappa,Bmat_TOT);
    end
    % correct nearfield
    IND_attached=find(tA_sub);
    S=[S_1m1_fast(1:3*N,:),S_1m1_fast(3*N+1:6*N,:),S_1m1_fast(6*N+1:9*N,:)];
    
    for ind=IND_attached'
        i=mod(ind-1,n_sig)+1;
        j=(ind-i)/n_sig+1;
        sig_dense.indices=I(i);
        tau_dense.indices=J(j);
        S([i,i+n_sig,i+2*n_sig],[j,j+n_tau,j+2*n_tau])=SLP_pld_SLOW(H,sig_dense,tau_dense,assembler,kappa);
    end
    A=S;
    B=1;
    mat.A=A;
    mat.B=B;
    mat.LR=0;
    mat.mem=16.*size(A,1)*size(A,2);
elseif(strcmp(assembler.method,'aca'))
    col=@(j)BF_clust_one_dof(XI,ETA,g_sig,g_tau,B,j,kappa);
    row=@(i)BF_clust_one_dof(ETA,XI,g_tau,g_sig,B',i,kappa);
    [A,B,flag]=aca(row,col,3*n_sig,3*n_tau,assembler.tol);
    if flag
        assembler_flag.method='dense';
        assembler_flag.QUAD=assembler.QUAD;
        assembler_flag.B=assembler.B;
        mat=SLP_pld(H,block,assembler_flag,kappa);
        mat.LR=0;
    else
    mat.A=A;
    mat.B=B;
    mat.LR=1;
    mat.mem=16.*(size(A,1)+size(B,1))*size(A,2);
    end
else
    error('assembly method unknown');
end

end


function I=BF_FAR(XI,eta,g_sig,g_tau,kappa,Bmat)
    N=size(XI,1)/3;
    r=zeros(N,size(XI,2));
    for i = 1:N
        r(i,:) = vecnorm(eta-XI(3*(i-1)+1:3*i,:));
    end    
    E=exp(1i*kappa*r);
    I=(E.*((g_tau.*g_sig)./(4*pi*r)))*Bmat;
    I=I(:);%reshape(I,[9*N,1]);
end


function I=BF_clust_one_dof(XI,ETA,g_sig,g_tau,B,ind_col,kappa)
M=size(ETA,1)/3;

ind_eta=mod(ind_col-1,M)+1;
eta=ETA(3*(ind_eta-1)+1:3*ind_eta,:);
g_tau=g_tau(ind_eta);

dof=(ind_col-ind_eta)/M + 1;
Bmat=[B{1,dof},B{2,dof},B{3,dof}];
I=BF_FAR(XI,eta,g_sig,g_tau,kappa,Bmat);
end
