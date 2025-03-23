function [Ker] = AniGau(Scale,theta,xscaling,yscaling)
%ANIGAU 此处显示有关此函数的摘要
%   Generating a Scale×Scale anisotropic Gaussian kernel, theta is the
%   spinning angle, xscaling and yscaling determing the shape of the
%   eclipse.
Lambda=diag([xscaling,yscaling])*[cos(theta) -sin(theta);sin(theta) cos(theta)];
Lambda=Lambda'*Lambda;
xy=[-floor((Scale-1)/2):floor(Scale/2)];
xy=[kron(xy,ones(1,Scale));kron(ones(1,Scale),xy)]';
Ker=sqrt(det(Lambda))*exp(reshape(diag(-xy*Lambda*xy'/2),Scale,[]))/2/pi;
% % % imshow(Ker,[])
end

