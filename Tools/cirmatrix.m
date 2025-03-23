function [cmatrix] = cirmatrix(v,mode,stride)
%CIRMATRIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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
