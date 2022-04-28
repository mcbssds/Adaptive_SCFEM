function M = stochcol_getalphaset(n)
%STOCHCOL_GETALPHASET generates the set {0,1}^n
%
%   TIFISS function: FX 16 April 2019
% Copyright (c) 2019 F. Xu

M = [0;1];
vec = [0;1];
nn = 1;
while nn<n
M = stoch_matrix_append(M,vec);
nn = nn +1;
end
M = sortrows(M);