function L = stochcol_getinterpolant_2(grid, clincombiset, ...
    indlincombiset, yy, polys)
%STOCHCOL_GETINTERPOLANT_2 evaluates multivariate Lagrange polynomials at collocation nodes
%
% L = stochcol_getinterpolant_2(grid, clincombiset, indlincombiset, yy, polys)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

[gr, gc] = size(grid);
[yr, yc] = size(yy);
if gc ~= yc
    error('dimension of test point is not same as dimension of collocation point!')
end
L = zeros(gr, yr);
for n = 1:yr
for i = 1:gr
    L(i,n) = 0;
    pts = grid(i,:);
    for j = 1:length(clincombiset{i})
        temp = 1;
        index = indlincombiset{i}(j,:);
        for k = 1:gc
            temp = temp*polys{index(k)}{pts(k)}(yy(n,k));
        end
        L(i,n) = L(i,n) + clincombiset{i}(j)*temp;
    end
end
end
