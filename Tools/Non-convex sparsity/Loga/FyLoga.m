function [z] = FyLoga(y,x,lambda,gamma)
%FY 此处显示有关此函数的摘要
%   此处显示详细说明
z=1/2*(y-x).^2+lambda*log(gamma*x+1)./log(gamma+1);
end

