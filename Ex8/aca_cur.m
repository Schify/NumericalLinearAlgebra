%% ACA-CUR script
clear;clc;
n = 800;
m = 800;
k=600;
A0 = randn(n,k);
B0 = randn(m,k);
U= orth(A0);
V = orth(B0);
l = linspace(0,16,k);
s = 10.^(-l);
S=diag(s);
E = randn(n,m);
E = E/norm(E);
M = U*S*V';%+(1e-14)*E;
M=(M/norm(M));
errS = zeros(15,1);
errU = zeros(15,1);
errI = zeros(15,1);
errLR = zeros(15,1);
for t=1:15
    tol = 10^(-t);
    [A,B,Ipiv,Jpiv,flag]=aca(M,n,m,tol);    %todo by you
    [A0,B0,UI0,d0,UJ0,Ipiv0,Jpiv0,flag0]=aca_cur_stable(M,n,m,tol); %todo by you
    [A1,B1,UI1,d1,UJ1,Ipiv1,Jpiv1,flag1]=aca_cur_naive(M,n,m,tol);
    
    M00U = (M(:,Jpiv0)*UJ0*diag(d0)*UI0.'*M(Ipiv0,:));
    M01U = (M(:,Jpiv1)*UJ1*diag(d1)*UI1.'*M(Ipiv1,:));
    M0i = M(:,Jpiv)*(M(Ipiv,Jpiv)\M(Ipiv,:));
    M0lr=A*B.';
    
    errS(t) = norm(M00U-M);
    errU(t) = norm(M01U-M);
    errI(t) = norm(M0i-M);
    errLR(t) = norm(M0lr-M);
    disp(t)
end