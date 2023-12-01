function y=matvec_HMAT(HMAT,x,Hx,Hy)

x_hat=x;%Hx.T'*x;
N=Hx.grid.N;
y=zeros(3*N,1);
BLOCKS=HMAT.BLOCKS;
for b=1:size(BLOCKS,2)
    block=BLOCKS{b};

    tau=block.tau.indices;
    sig=block.sig.indices;
    if(block.adm)
    y([sig;sig+N;sig+2*N])=y([sig;sig+N;sig+2*N])+block.mat.A*(block.mat.B.'*x_hat([tau;tau+N;tau+2*N]));
    y([tau;tau+N;tau+2*N])=y([tau;tau+N;tau+2*N])+block.mat.B*(block.mat.A.'*x_hat([sig;sig+N;sig+2*N]));
    else
    y([sig;sig+N;sig+2*N])=y([sig;sig+N;sig+2*N])+block.mat.A*(block.mat.B.'*x_hat([tau;tau+N;tau+2*N]));
    end
end
end