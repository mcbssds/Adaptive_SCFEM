function Y = stochcol_getgrid(X)
% STOCHCOL_GETGRID compute the set of indices for all sparse grid terms
%   grid = stochcol_getgrid(X)
%
% This function also compute the minimal downward closed set including X.
%   TIFISS function: FX 10 October 2019
% Copyright (c) 2019 F. Xu

% prepare X
X = sortrows(X);
Y = stochcol_gettensorgrid(X(end,:));
if size(X,1) > 1
    for i = 2:size(X, 1)
        index = X(end - i + 1,:);
        if ~ismember(index, Y, 'rows')
            Y = unique([Y; stochcol_gettensorgrid(index)],'rows');
        end
    end
end
Y = sortrows(Y);