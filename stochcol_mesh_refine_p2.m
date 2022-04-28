function paras_fem_new = stochcol_mesh_refine_p2(MMele, MMedge, paras_fem, ...
    paras_detail, pmethod)
%STOCHCOL_MESH_REFINE_P2   error estimation for reference solution
xy    = paras_fem{1};
evt   = paras_fem{2};
bound = paras_fem{3};
evtY    = paras_detail{4};
xyY     = paras_detail{5};
boundY  = paras_detail{6};
if pmethod == 1
    [evt, xy, bound, ~, eboundt] = mesh_ref(MMele, MMedge, evt, xy, ...
        bound, evtY, xyY, boundY);
    paras_fem_new = {xy, evt, bound, eboundt};
    xyp1 = xy;
    evtp1 = evt;
    boundp1 = bound;
    [xy, evt, bound] = p2_grid_generator(xy, evt, bound);
    paras_fem_new = {xy, evt, bound, eboundt, xyp1, evtp1, boundp1};
elseif pmethod == 2
    xyp1 = paras_fem{5};
    evtp1 = paras_fem{6};
    boundp1 = paras_fem{7};
    % Fixing data if P2-approximations have been used
    xy    = xyp1;
    evt   = evtp1;
    bound = boundp1;
    [evt, xy, bound, ~, eboundt] = mesh_ref(MMele, MMedge, evt, xy, ...
        bound, evtY, xyY, boundY);
    % Save oldest xy and bound data
    xyp1 = xy;
    evtp1 = evt;
    boundp1 = bound;
    [xy, evt, bound] = p2_grid_generator(xy, evt, bound);
    paras_fem_new = {xy, evt, bound, eboundt, xyp1, evtp1, boundp1};
end
