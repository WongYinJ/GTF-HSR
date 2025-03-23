function [G,Gr,ob] = G_Coding(X,Y,P,G,U,lambda,rho,nu,turank,shape_re,MaxLoop,ConvTol)
%RETF 此处显示有关此函数的摘要
%   此处显示详细说明
% X:HSI; Y: MSI; P:R×3 cell containing degradation matrices
% eta: low-rank parameter; phi: sparse parameter; rho: optimization penalty
%%%Preparation
[m,n,S]=size(X);[M,N,s]=size(Y);
ratio=M/m;
R=size(P,1);
K=length(shape_re);

% % Gu=tucker_als(tensor(imresize(X,ratio)),turank);
% % G=double(Gu.core);
% % U=Gu.U';Uaux=Gu.U';
% % psnr(double(ttm(tensor(G),U)),SRI)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % U=cell(1,3);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % [U{1},~,~]=svds(Unfold(Y,1)*Unfold(Y,1)',turank(1));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % [U{2},~,~]=svds(Unfold(Y,2)*Unfold(Y,2)',turank(2));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Gu=tucker_als(tensor(imresize(X,ratio)),turank);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % [U{3},~,~]=svds(Unfold(X,3)*Unfold(X,3)',turank(3));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % U{3}=Gu.U{3};
% % psnr(double(ttm(tensor(G),U)),SRI)

Gr=cell(1,R);
Pr=cell(1,R);
for r=1:R
    Gr{r}=G;Pr{r}=zeros(size(G));
end
Px=zeros(size(X));Py=zeros(size(Y));Ps=zeros(size(G));%Ps=zeros(prod(shape_re(1:2:N)),prod(shape_re(2:2:N)));

rhomax=1e20;
ps=[];
disp('Main Loop Starting...')
%% Main algorithm
for loop=1:MaxLoop
    Gpre=G;
    
    
    %%%Update Gr
    for r=1:R
        Tau=X+Px/rho-FormHSI(Gr([1:r-1,r+1:R]),U,U,P([1:r-1,r+1:R],:));
        Gr{r}=Update_G(Tau,{P{r,1}*U{1},P{r,2}*U{2},U{3}},G+Pr{r}/rho,1);
    end
    %%%Update Sg
    
    
% % % % %
% Sg=Lap(BlockUnfold(G-Ps/rho,shape_re),2*lambda/rho,1e-7,1e-7);Former Version for pami paper



        Sg=LapSparse(BlockUnfold(G-Ps/rho,shape_re),2*lambda/rho,1e0,1e-7);
    %%%Update G
    Wg=0;
    for r=1:R
        Wg=Wg+Gr{r}-Pr{r}/rho;
    end
    Wg=Wg+BlockFold(Sg,shape_re)+Ps/rho;
    G=Update_G(Y+Py/rho,{U{1},U{2},P{1,3}*U{3}},Wg/(R+1),1);
    %%%Update Lagragian 
    Px=Px+rho*(X-FormHSI(Gr,U,U,P));
    Py=Py+rho*(Y-double(ttm(tensor(G),{U{1},U{2},P{1,3}*U{3}})));
    for r=1:R
        Pr{r}=Pr{r}+rho*(G-Gr{r});
    end
    Ps=Ps+rho*(BlockFold(Sg,shape_re)-G);
    %%Update Penalty
    rho=min(nu*rho,rhomax);
    
    renewal=norm(Unfold(G-Gpre,1))/norm(Unfold(Gpre,1));
    Z=double(ttm(tensor(G),U));
% %     psnr(Z,SRI)
% %     ps=[ps,psnr(Z,SRI)];
    if mod(loop, 20) == 0 | loop==1
        fprintf('Iteration No.%d，Convergence criterion=%f\n', loop, renewal);
    end
    if renewal<ConvTol
        fprintf('======================Iteration No.%d，Algorithm converges======================', loop);
        break
    end
end
end

