function [X] = Fold(X,scale,dim)
%FOLD 此处显示有关此函数的摘要
%   此处
N=length(scale);
X=permute(reshape(X,[scale(dim),scale(1:dim-1),scale(dim+1:N)]),[2:dim,1,dim+1:N]);
end

