function [Unfolded] = BlockUnfold(G,shape_re)
%BLOCKUNFOLD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N=length(shape_re);
Unfolded=reshape(permute(reshape(G,shape_re),[1:2:N,2:2:N]),prod(shape_re(1:2:N)),[]);
end

