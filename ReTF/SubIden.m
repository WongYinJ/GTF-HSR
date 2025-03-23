function [U] = SubIden(X,Pc,Uy,s2,mu,rho,Ac,IniLoop)
%Spatial subspace identification
% Bc=Ac;
R=size(Ac,1);
A=[];
Pr=[];
for r=1:R
    A=[A;Ac{r}];
    Pr=[Pr,Pc{r}];
end
B=A;
M=zeros(size(A));
Ux=zeros(size(Uy,1),s2);
s1=size(Uy,2);

rhomax=1e20;
nu=1.3;
for loop=1:IniLoop
    Uxpre=Ux;
    Apre=A;
    Bpre=B;
    Mpre=M;
    
    % Update Ux
    Xy=X;
    for r=1:R
        Bry=B((r-1)*(s1+s2)+1:(r-1)*(s1+s2)+s1,:);
        Xy=Xy-Pc{r}*Uy*Bry;
    end
    E=cell(1,R^2);
    F=E;
    k=1;
    for r1=1:R
        for r2=1:R
            E{k}=Pc{r1}'*Pc{r2};
            F{k}=B(r2*(s1+s2)-s2+1:r2*(s1+s2),:)*B(r1*(s1+s2)-s2+1:r1*(s1+s2),:)';
            k=k+1;
        end
    end
    Res=0;
    for r=1:R
        Res=Res+Pc{r}'*Xy*B(r*(s1+s2)-s2+1:r*(s1+s2),:)';
    end
    Ux = CGra(E,F,ones(1,R^2),Ux,Res);
    %Update B
    D=sparse(Pr)*kron(sparse(eye(R)),sparse([Uy,Ux]));
    B=(D'*D+rho/2*eye(size(D,2)))^(-1)*(D'*X+rho/2*(A+M/rho));
    % Update A
    A=soft(B-M/rho,mu/rho);
    % Update M
    M=M+rho*(A-B);
    
    rho=min(nu*rho,rhomax);
    
    a=norm(Ux-Uxpre,'fro')
    b=norm(A-Apre,'fro')
    c=norm(B-Bpre,'fro')
    d=norm(M-Mpre,'fro')
    if max([a,b,c,d])<2e-4
        break
    end
end
U=[Uy,Ux];
end

