function A=HMAT(H_pld,bct,kappa,assembler)
% NxM H-matrix A
LEAVES=bct.leaves;


nL=size(LEAVES,2);
nB=0;
for i = 1:nL
    if(min(LEAVES{i}.sig.indices)<=min(LEAVES{i}.tau.indices))
     nB=nB+1;
    end
 end
BLOCKS=cell(1,nB);
mem_adm=0;
mem_dense=0;
block_ind=0;
for b=1:nL
    block=LEAVES{b};
    if(block.adm && min(block.sig.indices)<min(block.tau.indices))
        block_ind=block_ind+1;
        BLOCKS{block_ind}.sig=block.sig;
        BLOCKS{block_ind}.tau=block.tau;
        mat=SLP_pld(H_pld,block,assembler,kappa);
        BLOCKS{block_ind}.mat=mat;
        BLOCKS{block_ind}.LR=mat.LR;
        BLOCKS{block_ind}.adm=1;
        mem_adm=mem_adm+BLOCKS{block_ind}.mat.mem;
    elseif(~block.adm)
        block_ind=block_ind+1;
        block=LEAVES{b};
        BLOCKS{block_ind}.sig=block.sig;
        BLOCKS{block_ind}.tau=block.tau;
        mat=SLP_pld(H_pld,block,assembler,kappa);
        BLOCKS{block_ind}.mat=mat;
        BLOCKS{block_ind}.LR=0;
        BLOCKS{block_ind}.adm=0;
        mem_dense=mem_dense+BLOCKS{block_ind}.mat.mem;
    end
    if mod(b,100)==0
    blockstring=strcat(num2str(nL)," Blocks, ");
    str=strcat(strcat(strcat(strcat( num2str(100.*b/nL,'%.2f'),"% of "),blockstring) , strcat(num2str(mem_adm/1e9,'%.2f'),' GB ,') ),strcat(num2str(mem_dense/1e9,'%.2f'),' GB'));
    %str=strcat(strcat(num2str(b),'/'),num2str(nL));
    disp(str)
    end
end
A.mem=(mem_dense+mem_adm)/1e9; %mem in GB
A.BLOCKS=BLOCKS;
end