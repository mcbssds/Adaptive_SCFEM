function [nodes] = stochcol_nodes_cc(levels)
%STOCHCOL_NODES_CC computes the Clenshaw-Curtis nodes in [-1,1]
%
% input:
%     levels    vector of levels
%
% output:
%     nodes     cell of Clenshaw-Curtis collocation points
%
%   TIFISS function: FX 15 April 2019
% Copyright (c) 2019 F. Xu

nodes = cell(length(levels),1);
n = max(levels);
[nodes_increasing] = stochcol_nodes_cc_increasing(1:n);
% generate cc points sequence
xcc = nodes_increasing{1};
for k = 2:n
    xcc = [xcc setdiff(nodes_increasing{k},nodes_increasing{k-1})];
end
% choose set of cc points by level
for i=1:length(levels)
    if levels(i) == 1
        nodes{i} = xcc(1);
    else
        nodes{i} = xcc(1:(2^(levels(i)-1)+1));
    end
end
end

