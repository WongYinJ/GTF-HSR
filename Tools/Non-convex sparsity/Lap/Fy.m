function [z] = Fy(y,x,lambda,gamma)
%FY �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
z=1/2*(y-x).^2+lambda*(1-exp(-x/gamma));
end

