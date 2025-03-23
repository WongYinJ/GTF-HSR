function [Unfolded] = BlockFold(G,Folded_shape)
%BLOCKUNFOLD �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N=length(Folded_shape);
perdim=zeros(1,N);
perdim(1:2:N)=1:N/2;perdim(2:2:N)=N/2+1:N;
Unfolded=reshape(permute(reshape(G,[Folded_shape(1:2:N),Folded_shape(2:2:N)]),perdim),Folded_shape(1:2:N).*Folded_shape(2:2:N));
end

