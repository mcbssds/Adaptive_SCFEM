function [nodes] = stochcol_nodes_leja(levels)
%STOCHCOL_NODES_LEJA computes the Leja nodes in [-1,1]
%
% input:
%     levels    vector of levels
%
% output:
%     nodes     cell of Leja collocation points
%
% The function returns the Leja interpolation sparse grid nodes for a given
% set of levels 'i' (i>0).
%
% At level i, the number of collocation points are:
%
%   mi := i
%
% Given the first Leja point X(1) = 1, the Leja points are defined
% recursively in such a way that the function
%
% F(Y,X(1:N-1)) = ABS(PROD(Y-X(1:N-1)))
%
% is maximized over [-1,1] by Y = X(N). When more than one Y give the
% same maximum of the function, we choose the smallest as X(N). The maximum
% of the function aboveis computed by FMINBND (the algorithm is based on 
% golden section search and parabolic interpolation).
%
% This function is based on lejapoints.m by Marco Caliari.
% Available at http://profs.scienze.univr.it/~caliari/software.htm
%
%   TIFISS function: FX 15 April 2019
% Copyright (c) 2019 F. Xu

nodes = cell(length(levels),1);
n = max(levels);
options = [];
xleja = zeros(n,1);
xmaxloc = zeros(n-2,1);
fmaxloc = xmaxloc;
% generate leja points sequence
xleja(1) = 1;
xleja(2) = -1;
xsort = xleja;
for k = 3:n
    xsort = sort(xsort(1:k-1));
    for i = 1:(k-2)
        [xmaxloc(i),fmaxloc(i)] = fminbnd(@fleja,xsort(i),xsort(i+1),options,xleja,k-1);
    end
    [fmax index] = max(-fmaxloc(1:k-2));
    xleja(k) = xmaxloc(index);
    xsort(k) = xleja(k);
end
xleja = xleja(1:n);
% choose set of leja points by level
for i = 1:length(levels)
    nodes{i} = xleja(1:levels(i))';
end


    function z = fleja(x,xleja,n)
        % the opposite of the Leja function to maximize
        z = -abs(prod(x-xleja(1:n)));
    end

    function [min minval] = fminbnd(func,lb,ub, options, varargin)
        delta = 1e-17;
        gr = (sqrt(5)-1)/2;
        width = (ub-lb);
        out = [ lb:(width/3):ub ];
        out(2) = out(4)-gr*width;
        out(3) = out(1)+gr*width;
        upper = feval(func,out(3), varargin{:});
        lower = feval(func,out(2), varargin{:});
        while((out(3)-out(2)) > delta)
            if (upper > lower)
                out(4) = out(3);
                out(3) = out(2);
                width = out(4)-out(1);
                out(2) = out(4)-gr*width;
                upper = lower;
                lower = feval(func,out(2), varargin{:});
            else
                out(1) = out(2);
                out(2) = out(3);
                width = out(4)-out(1);
                out(3) = out(1)+width*gr;
                lower = upper;
                upper = feval(func,out(3), varargin{:});
            end
        end
        min = out(2);
        minval = feval(func,out(2), varargin{:});
    end
end