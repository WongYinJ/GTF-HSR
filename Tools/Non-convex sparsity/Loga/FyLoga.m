function [z] = FyLoga(y,x,lambda,gamma)
%FY �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
z=1/2*(y-x).^2+lambda*log(gamma*x+1)./log(gamma+1);
end

