function [lagpoly, weight]= stochcol_1Dlagpoly(node,nodes)
%STOCHCOL_1DLAGPOLY function handle and weight of 1D Lagrange polynomial
%
% [lagpoly, weight]= stochcol_1Dlagpoly(node,nodes)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

if length(nodes)==1
    lagpoly = @(y) 1;
    weight = 2;
else
    nodes = setdiff(nodes,node);
    n = length(nodes);
    lagpoly = @(y) 1;
    for k = 1:n
        lagpoly = @(y) lagpoly(y).*(y-nodes(k))./(node-nodes(k));
    end
    weight = integral(lagpoly,-1,1);
end
end
