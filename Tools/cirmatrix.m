function [cmatrix] = cirmatrix(v,mode,stride)
%CIRMATRIX 此处显示有关此函数的摘要
%   此处显示详细说明
i=0;
cmatrix=[];
if mode
    while i<length(v)
        cmatrix=[cmatrix;circshift(v,i)];
        i=i+stride;
    end
else
    while i<length(v)
        cmatrix=[cmatrix,circshift(v,i)];
        i=i+stride;
    end 
end
