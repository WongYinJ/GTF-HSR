function [z] = J2(y,x,lambda,gamma)
%J2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
t=J1(y,x,lambda,gamma);
tt=J1(y,t,lambda,gamma);
z=t-(tt-t).*(t-x)/(tt-2*t+t);
end

