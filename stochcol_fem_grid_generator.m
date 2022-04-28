function paras_fem = stochcol_fem_grid_generator(pmethod, dom_paras)
%STOCHCOL_FEM_GRID_GENERATOR triangular mesh generator
%
%  paras_fem = stochcol_fem_grid_generator(pmethod, dom_paras)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

dom_type = dom_paras(1);
if dom_type == 1 % square domain
    [~, ~, xy, mv, bound, mbound] = square_domain_fa(1, 1); % Q2 grid
    [evt, eboundt] = p1grid(xy,mv,bound,mbound,0); % P1 grid based on Q2 grid
elseif dom_type == 2 % L-shaped domain
    mesh_type = dom_paras(2);
    if mesh_type == 1
        [~, ~, xy, mv, bound, mbound] = ell_domain_fa;
        [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
    elseif mesh_type == 2
        [~, ~, xy, evt, bound, eboundt] = ell_domain_unstructured_fa;
        load ell_grid;
        % For unstructured meshes reorder the nodes of each element in evt
        [evt, xy, eboundt] = adjust_unstruct_mesh(evt,xy,eboundt);
    else
        error('Invalid mesh type');
    end
else
    error('Invalid mesh type');
end
if pmethod == 1
    paras_fem = {xy, evt, bound, eboundt};
elseif pmethod == 2
    % Save oldest xy and bound data
    xyp1 = xy;
    evtp1 = evt;
    boundp1 = bound;
    [xy, evt, bound] = p2_grid_generator(xy, evt, bound); % P2 grid
    paras_fem = {xy, evt, bound, eboundt, xyp1, evtp1, boundp1};
end
