function [Unfolded] = BlockUnfold(G,shape_re)
%BLOCKUNFOLD 此处显示有关此函数的摘要
%   此处显示详细说明
N=length(shape_re);
Unfolded=reshape(permute(reshape(G,shape_re),[1:2:N,2:2:N]),prod(shape_re(1:2:N)),[]);
end

