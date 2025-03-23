clear all;clc;
addpath(genpath('.\Urban'));
addpath(genpath('.\Tools'));
addpath(genpath('.\ReTF'));


SRI = multibandread('URBAN', [307,307,210], 'int16', 0, 'bil', 'ieee-be'); % SRI reference
HyWave=importdata('HydiceWave.mat');
SRI=SRI(1:256,1:256,:);
SRI(:,:,[1:4,76,87,101:111,136:153,198:210])=[];
HyWave([1:4,76,87,101:111,136:153,198:210],:)=[];
SRI=SRI.^(0.8);
SRI = SRI./max(SRI(:));

[I, J, K] = size(SRI);
S1=reshape(SRI,[I,K*J])'; % mode 1 matricization
S3=reshape(SRI,[I*J,K]); % mode 3 matricization
% ============================================
%% --------generate MS image
Landwavelength = [450 520; 520 600; 630 690; 760 900; 1550 1750; 2080 2350]; % ETM/Landsat
K_M=size(Landwavelength,1);
PM = zeros(K_M,K);
for b=1:K_M;
    b_i = find(HyWave(:,2)>Landwavelength(b,1),1);
    b_e = find(HyWave(:,2)<Landwavelength(b,2),1,'last');
    PM(b,b_i:b_e) = 1/(b_e+1-b_i);
end
% PM=sparse(PM);


M3=S3*PM';
MSI = reshape(M3,I,J,[]); % MSI
% K_M=size(MSI,3);

%% ===========================================
%constract HS image
ratio = 8; % downsampling ratio
Ih=I/ratio;
Jh=J/ratio;
kernel_length=ratio+1;
Ker=AniGau(kernel_length,3*pi/16,0.3,1.2);%% Anisotropic Gaussian kernel
s0=1;kertol=1;
[HSI,Pr,D] = SpatialDegrad(SRI,Ker,ratio,s0,kertol,'center'); % create HSI
% H3=reshape(HSI,[Ih*Jh,K]);
%%
disp('======================BGS-GTF-HSR applied to Anisotropic Gaussian blurring======================')
P=Pr;P{1,3}=PM;
rho=1e-1;nu=1.1;lambda=1e4;mu=3e0;
turank=[256,256,12];addatoms=[0,0];shape_re=[8,32,8,32,3,4];
MaxLoop=500;IniLoop=1000;
ConvTol=1e-5;%% Convergence Tolerance
Z_ReTF=ReTF(HSI,MSI,P,lambda,mu,rho,nu,turank,addatoms,shape_re,IniLoop,MaxLoop,ConvTol);
[psnr_LRTA,rmse_LRTA, ergas_LRTA, sam_LRTA, uiqi_LRTA,ssim_LRTA,DD_LRTA,CC_LRTA] = quality_assessments(SRI*255, Z_ReTF*255, 0, 1.0/ratio);
metrics=[psnr_LRTA,rmse_LRTA, ergas_LRTA, sam_LRTA, uiqi_LRTA,ssim_LRTA,DD_LRTA,CC_LRTA]