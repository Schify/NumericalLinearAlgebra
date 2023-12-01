%% Hmat exercise

load("HmatData.mat");
P=grid.P;
tri=grid.tri;
dir=[1;1;1];
dir=dir/norm(dir);
nvec=grid.n*dir;

figure(1)
clf
colormap jet
patch('Faces',tri,'Vertices',P,'FaceVertexCData',nvec,'FaceColor','flat')
bcw = readmatrix('bent-cool-warm-table-float-1024.csv');
bcw = bcw(:,2:4);
colormap(bcw)
axis equal
axis off
%%

kappa=.1*(grid.kappamax);
S_HMAT=HMAT(H_plc,bct,kappa,assembler_HMAT);
%%
d=[1.;1.;1.;];
d=d/norm(d);

u_inc=@(p)(exp(1i*kappa*(p*d)));
u_inc_vec=u_inc(P); %interpolatory
rhs=-u_inc_vec;
SLP_1m1=@(x) matvec_HMAT(S_HMAT,x,H_plc);%SLP_pld(H_pld,sig,tau,assembler,kappa);
T=H_plc.T;
SLP_1=@(x) T*SLP_1m1(T'*x);
Mass=H_plc.Mass;
MS=@(x)Mass\SLP_1(x);
conv_tol = 1e-5;
tic
[x,FLAG,~,~,res]=gmres(MS , rhs ,200, conv_tol , 1);
toc
%%
figure(2)
clf
patch('Faces',tri,'Vertices',P,'FaceVertexCData',(imag(x)),'FaceColor','interp','EdgeColor','none')
colormap(bcw)
axis off
axis equal