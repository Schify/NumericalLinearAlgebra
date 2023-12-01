
function S=SLP_pld_SLOW(H_plc,sig,tau,assembler,kappa)
% Function that computes the Single Layer Potential Boundary operator
% w.r.t. piecewise linear DIScontinuous basis functions.
% SLOW AND EXACT, MEANT FOR COMPARISON, TESTING AND ATTACHED TRIANGLES


grid=H_plc.grid;
P=grid.P;
tri=grid.tri;

QUAD = assembler.QUAD;


I=sig.indices';
J=tau.indices';

n_sig=numel(I);
n_tau=numel(J);
S=zeros(3*n_sig,3*n_tau);

for i_local=1:numel(I)
    i=I(i_local);
    p=tri(i,:)';
    for j_local=1:numel(J) 
        j=J(j_local);
        q=tri(j,:);
        for dofj=0:2
            S(i_local,j_local+dofj*n_tau)=BF_HH_pldc_new(p,q',P,QUAD,0,dofj,kappa,i,j);
            S(i_local+n_sig,j_local+dofj*n_tau)=BF_HH_pldc_new(p,q',P,QUAD,1,dofj,kappa,i,j);
            S(i_local+2*n_sig,j_local+dofj*n_tau)=BF_HH_pldc_new(p,q',P,QUAD,2,dofj,kappa,i,j);
        end            
    end
    %disp(i)
end
end