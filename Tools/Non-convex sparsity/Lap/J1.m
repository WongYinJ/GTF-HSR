function [z] = J1(y,x,lambda,gamma)
%J1 此处显示有关此函数的摘要
%   此处显示详细说明
z=y-lambda/gamma*exp(-x/gamma);
end

