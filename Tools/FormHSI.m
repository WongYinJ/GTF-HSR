function [HSI] = FormHSI(Gr,U,Uaux,P)
%FORMHSI �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
R=size(P,1);
HSI=0;
if length(Gr)>1
    for r=1:R
        HSI=HSI+double(ttm(tensor(Gr{r}),{P{r,1}*Uaux{1},P{r,2}*Uaux{2},U{3}}));
    end
else
    for r=1:R
        HSI=HSI+double(ttm(tensor(Gr{1}),{P{r,1}*Uaux{1},P{r,2}*Uaux{2},U{3}}));
    end
end
end

