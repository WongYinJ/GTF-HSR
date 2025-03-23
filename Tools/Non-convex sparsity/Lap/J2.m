function [z] = J2(y,x,lambda,gamma)
%J2 此处显示有关此函数的摘要
%   此处显示详细说明
t=J1(y,x,lambda,gamma);
tt=J1(y,t,lambda,gamma);
z=t-(tt-t).*(t-x)/(tt-2*t+t);
end

