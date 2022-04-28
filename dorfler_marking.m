function [MA, MB] = dorfler_marking(A, B, r)
%DORFLER_MARKING dorfler marking
%
% [MA, MB] = dorfler_marking(A, B, r)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

% sort contributions in descending order...
[B, sort_index] = sort(B,'descend');
% ... and sort index set accordingly
A = A(sort_index,:);
% sum of contributions
Bsum = sum(B);
% sum of contributions from marked indices 
MBsum = 0;
k = 1;
while MBsum < r*Bsum
    MBsum = sum(B(1:k));
    k = k + 1;
end
MB = B(1:k-1);
MA = A(1:k-1,:);