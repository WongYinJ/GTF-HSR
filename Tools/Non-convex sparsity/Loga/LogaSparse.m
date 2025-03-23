function [Smin] = LogaSparse(Y,lambda,gamma,tol)
%UNTITLED 此处显示有关此函数的摘要
%  To solve the problem of min 1/2*||S-Y||_F^2+lambda*\Sigma_i G(||s_i||2)
%  where G(x) is Logarithm
yaux=sqrt(sum(Y.^2,1));
index1=yaux~=0;
yau=yaux(:,index1);
%%%%
a0=max(-log((1-exp(-gamma))/lambda/gamma^2)/gamma,0);
index2=a0+lambda*gamma/(gamma*a0+1)/log(gamma+1)-yau<0;
y=yau(:,index2);
%%%%
% x=j*ones(size(y));
xi=y;
xiau=zeros(size(xi));
% yi=y;

iter=1;
maxIter=1e3;
while iter<maxIter
    indextol=abs(xi-xiau)>tol;
    xiau=xi;
%     indextol=abs(J1(y,J1(y,xi,lambda,gamma),lambda,gamma)-2*J1(y,xi,lambda,gamma)+xi)>tol;
%          abs(J1(y,J1(y,xi,lambda,gamma),lambda,gamma)-2*J1(y,xi,lambda,gamma)+xi)
% %     xiau=xi;xiau(~indextol)=j;x(x==j)=xiau;
% %     xi=xi(~indextol);yi=yi(~indextol);
    if sum(indextol)==0
        break
    end
    xi(indextol)=J1Loga(y(indextol),xi(indextol),lambda,gamma);
     iter=iter+1;
end
% abs(J1(y,J1(y,x,lambda,gamma),lambda,gamma)-2*J1(y,x,lambda,gamma)+x)
% % % % % % % % iter
x=J1Loga(y,xi,lambda,gamma);
xin2=zeros(size(yau));
xin2(index2)=x;xin2(~index2)=a0;
% indexFy=Fy(yau,xin2,lambda,gamma)>=Fy(yau,0,lambda,gamma);
xin2(FyLoga(yau,xin2,lambda,gamma)>=FyLoga(yau,0,lambda,gamma))=0;
Sin2=Y(:,index1)*sparse(diag(xin2./yau));
Smin=zeros(size(Y));
Smin(:,index1)=Sin2;
end

