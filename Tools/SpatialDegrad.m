function [I_m1,Pr,P] = SpatialDegrad(Im,Ker,ratio,s0,kertol,startpoint)
%SPATIALDEGRAD 此处显示有关此函数的摘要
%   This function blurrs image 'Im' using kernel 'Ker' and downsamples the
%   blurred result according to downsampling ratio 'ratio' and starting
%   point 's0', and returns the degraded image 'I_m', and the degradation
%   matrices in tensor factorization form 'Pr' and matrix form 'P'.
%   'kertol': to decide how many singular values of Ker are to be selected
[M,N,Dim]=size(Im);
FFTKer=psf2otf(Ker,[M,N]);
I_m=real(ifft2(fft2(Im).*Fold(repmat(FFTKer(:),[1,Dim])',[M,N,Dim],3)));
I_m=I_m(s0:ratio:end,s0:ratio:end,:);
%%% Kernel Truncation
[U,Sig,V]=svd(Ker);
sig=diag(Sig);
power=sig(1);
for R=1:min(size(Ker))
    if power/sum(sig)<kertol
        R=R+1;
        power=power+sig(R);
    else
        break
    end
end
%%%%%
Pr=cell(R,3);
P=sparse(0);
I_m1=0;
if strcmpi(startpoint,'center')
for r=1:R
    Pr{r,1}=cirmatrix(circshift(padarray(U(:,r)'*sqrt(sig(r)),[0 M-size(Ker,1)],'post'),-floor((size(Ker,1)+1)/2)+1),1,1);
%     Pr{r,1}(end-9:end,1:9)=0;Pr{r,1}(1:9,end-9:end)=0;
    Pr{r,1}=Pr{r,1}(s0:ratio:end,:);
    Pr{r,2}=cirmatrix(circshift(padarray(V(:,r)'*sqrt(sig(r)),[0 N-size(Ker,2)],'post'),-floor((size(Ker,2)+1)/2)+1),1,1);
%     Pr{r,2}(end-9:end,1:9)=0;Pr{r,2}(1:9,end-9:end)=0;
    Pr{r,2}=Pr{r,2}(s0:ratio:end,:);
    I_m1=I_m1+double(ttm(tensor(Im),{Pr{r,1},Pr{r,2}},[1 2]));
    P=P+kron(sparse(Pr{r,2}),sparse(Pr{r,1}))';
end
end
if strcmpi(startpoint,'first')
for r=1:R
    Pr{r,1}=cirmatrix(circshift(padarray(U(:,r)'*sqrt(sig(r)),[0 M-size(Ker,1)],'post'),0),1,1);
    Pr{r,1}=Pr{r,1}(s0:ratio:end,:);
    Pr{r,2}=cirmatrix(circshift(padarray(V(:,r)'*sqrt(sig(r)),[0 N-size(Ker,2)],'post'),0),1,1);
    Pr{r,2}=Pr{r,2}(s0:ratio:end,:);
    I_m1=I_m1+double(ttm(tensor(Im),{Pr{r,1},Pr{r,2}},[1 2]));
    P=P+kron(sparse(Pr{r,2}),sparse(Pr{r,1}))';
end
end
I_m2=Fold(Unfold(Im,3)*P,[M/ratio,N/ratio,Dim],3);
end

