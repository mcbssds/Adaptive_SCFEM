function G = stochcol_gmatrices(grid,clcset,indlcset,list)
%STOCHCOL_GMATRICES G-matrix calculator
%
%  G = stochcol_gmatrices(grid,clcset,indlcset,list)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

K = size(grid, 1);
N = size(grid, 2);
G = zeros(K, K);
for k = 1:K
    for l = 1:k
        temp = 0;
        cset_k = clcset{k};
        cset_l = clcset{l};
        indset_k = indlcset{k};
        indset_l = indlcset{l};
        for i = 1:length(cset_k)
            for j = 1:length(cset_l)
                temp_int = 1;
                for n = 1:N
                    temp_int = temp_int*list{indset_k(i,n), ...
                        indset_l(j,n)}(grid(k,n),grid(l,n));
                end
                temp = temp + temp_int*cset_k(i)*cset_l(j);
            end
        end
        G(k,l) = temp;
        G(l,k) = temp;
    end
end