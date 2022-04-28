function M = stoch_matrix_append(M1,vec)
% STOCH_MATRIX_APPEND add a vector after padding a matrix
% M = stoch_matrix_append(M1,vec);
%   SIFISS function: FX 8 March 2019
% Copyright (c) 2019 F. Xu

n = length(vec);
[s, t] = size(M1);
M = zeros(s*n,t+1);
if s == 0
    M = vec;
else
    for i = 1:n
        M(((i-1)*s+1):(i*s),:) = [M1 vec(i)*ones(s,1)];
    end
end