function [z] = J1(y,x,lambda,gamma)
%J1 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
z=y-lambda/gamma*exp(-x/gamma);
end

