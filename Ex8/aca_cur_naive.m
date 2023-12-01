function [A,B,UI,d,UJ,Ipiv,Jpiv,flag]=aca_cur_naive(M,n,m,tol)
%   LOW rank compression by adaptive cross approximation
%   INPUT:  row,col:    functions that generate row and col
%           n,m:        col,row size
%           tol:        preset relative tolerance (e.g. 1e-4)
%   OUTPUT: A,B:        LR pencil
%           flag:       true on failure

flag=0;
I=1:n;
J=1:m;
dim = randi(2,1);
kmax=floor(min(n,m)/2);
Ipiv=[];
Jpiv=[];

A=zeros(n,10);
B=zeros(m,10);

%if(dim==1)
i=randi(n,1);
b=M(I(i),:).';
[~,j]=max(abs(b));
a=M(:,j);
[~,i]=max(abs(a));
b=M(i,:).';
d=1/a(i);
a=a/a(i);

A(:,1)=a;
B(:,1)=b;
I=I(I~=i);
J=J(J~=j);
Ipiv=[Ipiv;i];
Jpiv=[Jpiv;j];
res0=norm(a)*norm(b);
res=Inf;
k=1;
UI=1;
UJ=1;
D=diag(d);
while res>res0*tol && k<n
    if(mod(k,20)==0)
        A=[A,zeros(n,10)];
        B=[B,zeros(m,10)];
    end
    dim = randi(2,1);
    %if(dim==1)
        i_ind=randi(n-k,1);
        i=I(i_ind);
        %b=M(i,:).'-B*(A(i,:).');
        b=M(i,:).'-(M(i,Jpiv)*(UJ*D*UI.')*M(Ipiv,:)).';
        [mx,j]=max(abs(b(J)));
        %a=M(:,J(j))-A*(B(J(j),:).');
        a=M(:,J(j))-M(:,Jpiv)*(UJ*D*UI.')*M(Ipiv,J(j));
        [mx,i]=max(abs(a(I)));
        %b=M(I(i),:).'-B*(A(I(i),:).');
        b=M(I(i),:).'-(M(I(i),Jpiv)*(UJ*D*UI.')*M(Ipiv,:)).';
        if mx>1e-15
            m00 = M(I(i),J(j))-M(I(i),Jpiv)*(UJ*D*UI.')*M(Ipiv,J(j));
            
            a=a/a(I(i));
            res=norm(a)*norm(b);
            
            A(:,k+1)=a;
            B(:,k+1)=b;

            %v=diag(d)*UI.'*M(Ipiv,J(j));
            %w=diag(d)*UJ.'*M(I(i),Jpiv).';
            
            UJ0 = [UJ,-1.*UJ*(diag(d)*(UI.'*M(Ipiv,J(j))));zeros(1,size(UJ,2)),1];
            UI = [UI,-1.*UI*D*UJ.'*M(I(i),Jpiv).';zeros(1,size(UI,2)),1];
            UJ=UJ0;
            d=[d,1/m00];
            D=diag(d);
            k=k+1;
        else
            res=0;
        end
    Ipiv=[Ipiv;I(i)];
    Jpiv=[Jpiv;J(j)];
    I=I(I~=I(i));
    J=J(J~=J(j));
    
    
end
A=A(:,1:k);
B=B(:,1:k);
Ipiv = Ipiv(1:k);
Jpiv = Jpiv(1:k);
UI=UI(1:k,1:k);
UJ=UJ(1:k,1:k);
d=d(1:k);


if k==kmax
    flag=1;
    disp('FLAG')
end
end