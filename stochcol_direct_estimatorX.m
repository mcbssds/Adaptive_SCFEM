function serrest = stochcol_direct_estimatorX(X, ...
    paras_sg, paras_fem, list, pmethod, rhs_fun, aa)
%STOCHCOL_DIRECT_ESTIMATORX computes spatial error estimate for SC
% designed by David Silvester in February 2021 updated December 2021
global Y1 Y2
G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
paras_detail = stochcol_mesh_detail(paras_fem);
MMele_full = (1:size(paras_fem{2}, 1))';      % all elements
MMedge_full = (1:size(paras_detail{7}, 1))';  % all edges
paras_fem_new = stochcol_mesh_refine(MMele_full, MMedge_full, ...
    paras_fem, paras_detail, pmethod);
a_unit = @(x1, x2) ones(size(x1));
[A_unit,~] = stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
    a_unit, rhs_fun);
[A_unit_new,~] = stochcol_fem_setup(paras_fem_new{1}, paras_fem_new{2}, ...
    a_unit, rhs_fun);
coords = paras_sg{9};
xy = paras_fem{1};
xy1 = xy(:,1);
xy2 = xy(:,2);
xy_new = paras_fem_new{1};
xy_new1 = xy_new(:,1);
xy_new2 = xy_new(:,2);
sols = zeros(length(xy), size(coords, 1));
sols_new = zeros(length(xy_new), size(coords, 1));
diff_sols = zeros(size(sols_new));
for k = 1:size(coords, 1)
% save parameter values as global variables
Y1=coords(k,1); Y2=coords(k,2); %disp(['check27'])
    [x_gal, ~] = stochcol_fem_solver(coords(k, :), paras_fem, ...
        aa, rhs_fun);
    [x_gal_new, ~] = stochcol_fem_solver(coords(k, :), paras_fem_new, ...
        aa, rhs_fun);
    x_gal_interp = griddata(xy1, xy2, x_gal, xy_new1, xy_new2);
    sols(:, k) = x_gal;
    sols_new(:, k) = x_gal_new;
    x_gal_diff = x_gal_new - x_gal_interp;
    diff_sols(:, k) = x_gal_diff;
end
serrest = sqrt(sum(dot(G, diff_sols' * A_unit_new * diff_sols)));
return
