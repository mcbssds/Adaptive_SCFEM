function Y = stochcol_margin_set(X)
%STOCHCOL_MARGIN_SET computen the margin set of the downward closed set
%  Y = stochcol_margin_set(X)
%
%   TIFISS function: FX 10 October 2019
% Copyright (c) 2019 F. Xu

% prepare X
X = sortrows(X);
% test whether X is downward closed;
Xdc = stochcol_getgrid(X);
if ~isequal(X, Xdc)
    error('The input is not a downward closed set!')
end

[r, c] = size(X);
UnitSet = eye(c);
XX = [];
for i = 1:r
    index = X(i,:);
    XX = unique([XX; ones(c, 1)*index + UnitSet], 'rows');
end

% return the margin set
Y = setdiff(XX, X, 'rows');
end