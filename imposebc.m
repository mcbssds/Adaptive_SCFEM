function [As,bs] = imposebc(A,b,intern,bound,bc)
% IMPOSEBC apply Dirichlet boundary conditions
%   [As,bs] = imposebc(A,b,intern,bound,bc)
%   input
%        A stiffness matrix without bc
%        b rhs without bc
%        intern interior nodes
%        bound  boundary nodes
%        bc     boundary values
 
% update the rhs
bs = b(intern);
bs = bs - A(intern,bound)*bc;
% update the stiffness matrix
As = A(intern,intern);

