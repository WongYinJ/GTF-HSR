function [z] = J1Loga(y,x,lambda,gamma)
%J1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
z=y-lambda*gamma./(gamma*x+1)./log(gamma+1);
end

