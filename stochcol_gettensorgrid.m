function M = stochcol_gettensorgrid(index)
%STOCHCOL_GETTENSORGRID generate tensor product grid
%
% M = stochcol_gettensorgrid(index)
%
% input
%       index   vector
% output 
%       M       matrix whose entries satisfy 1 <= M(j, i) < index(i)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

n = length(index);
M = [];
for i = 1:n
M = stoch_matrix_append(M,(1:index(i))');
end
