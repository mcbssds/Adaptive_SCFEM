function fun_pvalue = tgauss(s, t, xl_v, yl_v, fun)
% STOCHCOL_GAUSS evaluate function at physical points for triangle elements
%       coeff = stochcol_gauss_coeff(s,t,xl_v,yl_v,yy)
%       input
%           s      x coordinate of Gaussian point in refence element
%           t      y coordinate of Gaussian point in refence element
%           xl_v   x coordinate of vertex in physical element
%           yl_v   y coordinate of vertex in physical element
%           fun    function handle
% Copyright (c) 2019 Feng Xu
nel = length(xl_v(:,1)); % number of elements
zero_v = zeros(nel,1); 
x1 = zero_v; 
x2 = x1;
% vectorized linear shape functions
[phi_e,~,~] = vtshape(s, t);
% calculate physical coordinates of points in all physical elements which
% are corresponding to the point (s, t) in the reference element using
% isoparametric relation
for ivtx = 1:3 % size(phi_e, 2) == 3 
    x1 = x1 + phi_e(:, ivtx) .* xl_v(:, ivtx);
    x2 = x2 + phi_e(:, ivtx) .* yl_v(:, ivtx);
end
fun_pvalue = fun(x1, x2); % fun_pvalue has the size (nel, 1)
