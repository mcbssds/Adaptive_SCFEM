function [A, b] = stochcol_fem_setup(xy, mv, coeff_fun, rhs_fun)
%STOCHCOL_FEM_SETUP finite element matrix system generator
%
%   [A, b] = stochcol_fem_setup(xy, mv, diff_fun, rhs_fun)
%
%   input
%                xy  nodal coordinate vector
%                mv  element mapping matrix
%          p_method  approximation method (1 - linear; 2 - quadratic)
%         coeff_fun  function handle for the diffusion coefficient
%           rhs_fun  function handle for the right-hand side
%   output
%          A         stiffness matrix
%          b         rhs vector
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function imposebc.
[nel, dim] = size(mv); % % no. of elements and no. of nodes per element
if dim == 3
    p_method = 1;
elseif dim == 6
    p_method = 2;
end
%if p_method == 1
%    fprintf('\nsetting up P1 matrices... \n ')
%elseif p_method == 2
%    fprintf('\nsetting up P2 matrices... \n ')
%end
x = xy(:,1);
y = xy(:,2);
n = length(x); % number of nodes
% Gauss point integration rules (3/7/19/73 Gaussian points)
ngpt = 19;
[s,t,wt] = triangular_gausspoints(ngpt);
% inner loop over elements
% (xl_v(i,j),yl_v(i,j)) the coordinate of vertex j of element i
for ivtx = 1:3 % each triangular element has 3 vertices(P1 or P2)
    xl_v(:,ivtx) = x(mv(:,ivtx));
    yl_v(:,ivtx) = y(mv(:,ivtx));
end
% initialise element stiffness matrices
Ae = zeros(nel,dim,dim);
be = zeros(nel,dim);
% loop over Gauss points
for igpt = 1:ngpt
    sigpt = s(igpt);
    tigpt = t(igpt);
    wtigpt = wt(igpt);
    %  evaluate derivatives etc, of FE basis
    [jac,invjac,phi,dphidx,dphidy] = tderiv(sigpt,tigpt,xl_v,yl_v); % P1
    if p_method == 2
        [phi,dphidx,dphidy] = tqderiv(sigpt,tigpt,xl_v,yl_v); % P2
    end
    [coeff] = tgauss(sigpt, tigpt, xl_v, yl_v, coeff_fun);
    [rhs] = tgauss(sigpt,tigpt,xl_v,yl_v, rhs_fun);
    for j = 1:dim
        for i = 1:dim
            Ae(:,i,j) = Ae(:,i,j) + ...
                wtigpt*coeff(:).*dphidx(:,i).*dphidx(:,j).*invjac(:);
            Ae(:,i,j) = Ae(:,i,j) + ...
                wtigpt*coeff(:).*dphidy(:,i).*dphidy(:,j).*invjac(:);
        end
    end
    for i = 1:dim
        be(:,i) = be(:,i) + wtigpt*rhs(:).*phi(:,i).*jac(:);
    end
end
% assembly element contributions into global matrix/vector
A = assembly(Ae,mv,n,mv,n);
b = assembly(be,mv,n);
% A = sparse(n,n);
% b = zeros(n,1);
% for krow=1:dim
%     nrow=mv(:,krow);
%     for kcol=1:dim
%         ncol=mv(:,kcol);
%         A = A + sparse(nrow,ncol,Ae(:,krow,kcol),n,n);
%     end
%     for els=1:nel
%         b(nrow(els),1)=b(nrow(els),1) + be(els,krow);
%     end
% end
end
