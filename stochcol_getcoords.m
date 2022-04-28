function coords = stochcol_getcoords(grid, rule_id)
%STOCHCOL_GETCOORDS get coordinates of collocation nodes in parametric space
%from their grids
%
%  coords = stochcol_getcoords(grid,rule_id)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu


% coords is one-to-one with grid
if rule_id == 1
    levels = max(max(grid));
    nodes = stochcol_nodes_leja(levels);
elseif rule_id == 2
    if max(max(grid)) == 1
        levels = max(max(grid));
    else
        levels = log2(max(max(grid)-1))+1;
    end
    nodes = stochcol_nodes_cc(levels);
end
[rg, cg]=size(grid);
coords = zeros(rg,cg);
for i = 1:rg
    for j = 1:cg
        coords(i,j) = nodes{1}(grid(i,j));
    end
end
end

