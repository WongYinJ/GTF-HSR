function [z] = J1Loga(y,x,lambda,gamma)
%J1 此处显示有关此函数的摘要
%   此处显示详细说明
z=y-lambda*gamma./(gamma*x+1)./log(gamma+1);
end

