function [X] = Unfold(X,dim)
%UNFOLD 此处显示有关此函数的摘要
%   此处显示详细说明
scale=size(X);
N=length(scale);
X=reshape(permute(X,[dim,1:dim-1,dim+1:N]),scale(dim),[]);
end

