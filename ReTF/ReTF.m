function [Z,ob] = ReTF(X,Y,P,lambda,mu,rho,nu,turank,addatoms,shape_re,IniLoop,MaxLoop,ConvTol)
%RETF 此处显示有关此函数的摘要
%   此处显示详细说明
% X:HSI; Y: MSI; P:R×3 cell containing degradation matrices
% eta: low-rank parameter; phi: sparse parameter; rho: optimization penalty
%%%Preparation
disp('Initializing......')
[m,n,S]=size(X);[M,N,s]=size(Y);
ratio=M/m;
R=size(P,1);
K=length(shape_re);

% % Gu=tucker_als(tensor(imresize(X,ratio)),turank);
% % G=double(Gu.core);
% % U=Gu.U';Uaux=Gu.U';
% % psnr(double(ttm(tensor(G),U)),SRI)
U=cell(1,3);
[U{1},~,~]=svds(Unfold(Y,1)*Unfold(Y,1)',turank(1));
[U{2},~,~]=svds(Unfold(Y,2)*Unfold(Y,2)',turank(2));
% Gu=tucker_als(tensor(imresize(X,ratio)),turank);
% U{3}=Gu.U{3};
% [U{3},~]=vca(Unfold(X,3),turank(3))
[U{3},~,~]=svds(Unfold(X,3)*Unfold(X,3)',turank(3));

if prod(addatoms)~=0
Ac=cell(1,R);
for r=1:R
    Ac{r}=randn(turank(1)+addatoms(1),size(Unfold(X,1),2));
end
U{1}=SubIden(Unfold(X,1),P(:,1)',U{1},addatoms(1),mu,rho,Ac,IniLoop);
Ac=cell(1,R);
for r=1:R
    Ac{r}=randn(turank(2)+addatoms(2),size(Unfold(X,2),2));
end
U{2}=SubIden(Unfold(X,2),P(:,2)',U{2},addatoms(2),mu,rho,Ac,IniLoop);
end
corerank=turank;
corerank(1:2)=corerank(1:2)+addatoms;
G=G_Coding(X,Y,P,zeros(corerank),U,lambda,rho,nu,corerank,shape_re,MaxLoop,ConvTol);
% % % % % % % % % % % % % % % % % % % % % % % % % % G=rand(turank);
Z=double(ttm(tensor(G),U));
end

