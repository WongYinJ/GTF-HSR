function [G] = Update_G(T,U,W,lambda)
%UPDATE_G 此处显示有关此函数的摘要
%   此处显示详细说明
% This function deals with the problem of    
% min_G ||T-[G;U_1,U_2,U_3]||_F^2+\lambda||G-W||_F^2
V=cell(1,3);Sig=cell(1,3);
for i=1:3
    [V{i},Sig{i},~]=svd(U{i}'*U{i});
end
Gaux=double(ttm(tensor(T),{U{1}',U{2}',U{3}'}))+lambda*W;
Gaux=double(ttm(tensor(Gaux),{V{1}',V{2}',V{3}'}));
Gaux=Gaux./(double(ktensor({diag(Sig{1}),diag(Sig{2}),diag(Sig{3})}))+lambda);
G=double(ttm(tensor(Gaux),{V{1},V{2},V{3}}));
end

