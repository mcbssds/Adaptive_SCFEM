function [err_est, serrest, perrest] = stochcol_direct_estimator(X, ...
    paras_sg, paras_fem, list, rule_id, pmethod, rhs_fun, aa)
%------------------------------------------ legacy code
G = stochcol_gmatrices(paras_sg{4}, paras_sg{5}, paras_sg{6}, list);
%X_diff = stochcol_margin_set(X);
X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X)); %------ reduced
paras_sg_diff = stochcol_sg(X_diff, rule_id);
G_diff = stochcol_gmatrices(paras_sg_diff{4}, paras_sg_diff{5}, ...
    paras_sg_diff{6}, list);
paras_detail = stochcol_mesh_detail(paras_fem);
MMele_full = (1:size(paras_fem{2}, 1))';      % all elements
MMedge_full = (1:size(paras_detail{7}, 1))';  % all edges
paras_fem_new = stochcol_mesh_refine(MMele_full, MMedge_full, ...
    paras_fem, paras_detail, pmethod);
a_mean = @(x1, x2) ones(size(x1));
[A_mean,~] = stochcol_fem_setup(paras_fem{1}, paras_fem{2}, ...
    a_mean, rhs_fun);
[A_mean_new,~] = stochcol_fem_setup(paras_fem_new{1}, paras_fem_new{2}, ...
    a_mean, rhs_fun);
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
serrest = sqrt(sum(dot(G, diff_sols' * A_mean_new * diff_sols)));
coords_diff = paras_sg_diff{9};
sols_diff = zeros(length(paras_fem{1}), size(coords_diff, 1));
parfor k = 1:size(coords_diff, 1)
    [x_gal, ~] = stochcol_fem_solver(coords_diff(k, :), paras_fem, ...
        aa, rhs_fun);
    sols_diff(:, k) = x_gal;
end
perrest = sqrt(sum(dot(G_diff, sols_diff' * A_mean * sols_diff)));
% debug
X_diff,G_diff,coords_diff,sum(sols_diff), perrest,serrest
err_est = serrest + perrest;
