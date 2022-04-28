% driver for running single level adaptive collocation code
% designed by David Silvester in July 2021 updated 16-9-21
stochcol_adaptive_global_settings;
delta = default('set the error tolerance (default is 6e-3)',6e-3);
adaptmax = default('set the number of adaptive steps (default is 36)',36);
cpointsmax = 7;
tot_err_est = Inf;
tot_err_direct = Inf;

% preallocation of arrays
    dof = nan(1,adaptmax);
    err_p_iter = nan(1,adaptmax);
    err_s_iter = nan(1,adaptmax);
    error_iter = nan(1,adaptmax);
    
    err_p_d_iter = nan(1,adaptmax);
    err_s_d_iter = nan(1,adaptmax);
    error_d_iter = nan(1,adaptmax);

% preallocation of cell data
    sols_iter = cell(1,adaptmax);
    paras_sg_iter = cell(1,adaptmax);
    paras_fem_iter = cell(1,adaptmax);
    
iter = 0; glevel = 0;
startLoopTime = tic;
while tot_err_direct >= delta && iter <= adaptmax
    if iter == 0 % First iteration step
        % Initial index set is for a single collocation point
        X = stochcol_getindexset(0, M);
        % Several attributes of the general sparse grid interpolation given
        % the index set and the sequence of 1D collocation nodes, such as
        % the underlying grid, the coordinates of the underlying grid, the
        % multivariate Lagrange polynomials (represented by the coefficient
        % list of the single-term multivariate Lagrange polynomials and the
        % set of indexing the degree of the single-term polynomials), etc.
        paras_sg = stochcol_sg(X, rule_id);
        gridd = paras_sg{4};
        clincombiset = paras_sg{5};
        indlincombiset = paras_sg{6};
        coords = paras_sg{9};
        figure(901)
        scatter(coords(:,1),coords(:,2),100,'o','filled');
        axis square
        box on
        grid on
        title('initial collocation nodes in first two directions')
        % Several attributes of the finite element: vertex coordinate
        % vector, element mapping matrix, boundary vertex vector, boundary
        % elment mapping matrix
        paras_fem0 = stochcol_fem_grid_generator(pmethod, dom_paras);
        paras_fem = paras_fem0;
        % Collection of edge information for flux jump computation and 
        % outputs of linear detail space Y grid-generator
        paras_detail0 = stochcol_mesh_detail(paras_fem0);
        paras_detail = paras_detail0;
        % plot the first mesh
        if pmethod == 1 % P1
            plot_mesh(paras_fem0{2}, paras_fem0{1}, 'initial FE mesh');
            pause(1)
        elseif pmethod == 2
            plot_mesh(paras_fem0{6}, paras_fem0{5}, 'initial FE mesh');
        end
        % We represent the general sparse grid interpolation of the finite
        % element approximation as a cell of vertex values associated with
        % collocation nodes in the order they are introduced
        sols = zeros(length(paras_fem0{1}), size(coords, 1));
        % For each collocation node, we calculate an estimate of the
        % energy error of the FE solution and a list of
        % marked elements and edges
        errests = nan(1,size(coords, 1));
        MMeles  = cell(1,size(coords, 1));
        MMedges  = cell(1,size(coords, 1));
        % use parfor to enable parallel computing
        for k = 1:size(coords, 1)
            % FE solution for a current collocation node
            [x_gal, ~] = stochcol_fem_solver(coords(k, :), paras_fem0, ...
                aa, rhs_fun);
            % Energy error of the FE solution and marked elements and edges
            [errest, ~, MMele, MMedge] = stochcol_fem_estimator(...
                pmethod, paras_fem_errest, paras_fem0, paras_detail0, ...
                x_gal, coords(k, :), aa, aax1, aax2, rhs_fun);
            sols(:, k) = x_gal;
            sols_ml{k} = x_gal;
            errests(k) = errest;
            MMeles{k} = MMele;
            MMedges{k} = MMedge;
        end
        % 2-norm of multivariate Lagrange polynomials
        L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
        % construct the detail index set ---- reduced margin
        X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X));
%        X_diff = stochcol_margin_set(X);
        % Attributes of the general sparse grid interpolation corresponding
        % to the new index set, which is the union of the original index 
        % set and the detail index set
        paras_sg_full = stochcol_sg([X; X_diff], rule_id);
        % set grid due to the inclusion of the detail index set
        % the grid is used to label the collocation nodes
        [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, ...
            paras_sg{4}, 'rows');
        coords_diff = paras_sg_full{9}(grid_diff_ind,:);
    else
        fprintf('\n\nIteration %i \n',iter)
        fprintf(['   spatial error indicator', ...
         ' is %10.4e \n'],serrest)
        fprintf(['parametric error indicator', ...
         ' is %10.4e \n'],perrest)
        if serrest < perrest % serrest_d < perrest_d
            fprintf('Parametric enrichment step ... new indices added \n');
            glevel = glevel+1;
            % Refined index set
            disp([Mset])
            X = stochcol_getgrid([X; Mset]);
            % Attributes of the general sparse grid interpolation based on
            % the refined index set
            paras_sg_new = stochcol_sg(X, rule_id);
            % Grid of collocation nodes based on previous index set
            gridd = paras_sg{4}; 
            % Grid and coordinates of collocation nodes based on the
            % refined index set
            gridd_new = paras_sg_new{4};
            coords_new = paras_sg_new{9};
            % FE solutions, energy error estimates, marked elements and
            % edges for collocation nodes  based on the refined index set
            % All FE solutions and part of estimates are from last iteration
            sols_new = zeros(length(paras_fem{1}), size(coords_new, 1));
            errests_new = zeros(1,size(coords_new, 1));
            MMeles_new  = cell(1,size(coords_new, 1));
            MMedges_new  = cell(1,size(coords_new, 1));
            [~, IA, IB] = intersect(gridd_new, gridd, 'rows');
            [~, IC, ID] = intersect(gridd_new, gridd_diff, 'rows');
            for k = 1:length(IA)
                % reuse the FE solutions and estimates from the last iteration
                sols_new(:, IA(k)) = sols(:, IB(k));
                errests_new(IA(k)) = errests(IB(k));
                MMeles_new{IA(k)} = MMeles{IB(k)};
                MMedges_new{IA(k)} = MMedges{IB(k)};
            end
            for k = 1:length(IC)
                % reuse the FE solutions from the last iteration
                sols_new(:, IC(k)) = sols_diff(:, ID(k));
            end
            % Compute new FE error estimates
            coords_temp = coords_diff(ID, :);
            sols_diff_temp = sols_diff(:, ID);
            errests_temp = zeros(1,length(ID));
            MMeles_temp = cell(1,length(ID));
            MMedges_temp = cell(1,length(ID));
            for k = 1:length(ID)
                x_gal = sols_diff_temp(:,k);
                [errest, ~, MMele, MMedge] = stochcol_fem_estimator(pmethod, ...
                    paras_fem_errest, paras_fem, paras_detail, x_gal, coords_temp(k, :), ...
                    aa, aax1, aax2, rhs_fun);
                errests_temp(k) = errest;
                MMeles_temp{k} = MMele;
                MMedges_temp{k} = MMedge;
            end
            for k = 1:length(IC)
                errests_new(IC(k)) = errests_temp(k);
                MMeles_new{IC(k)} = MMeles_temp{k};
                MMedges_new{IC(k)} = MMedges_temp{k};
            end
            paras_sg = paras_sg_new;
            gridd = paras_sg{4};
            clincombiset = paras_sg{5};
            indlincombiset = paras_sg{6};
            coords = paras_sg{9};
            sols = sols_new;
            errests = errests_new;
            MMeles = MMeles_new;
            MMedges = MMedges_new;
            % 2-norm of multivariate Lagrange polynomials
            L_two_norm = stochcol_multilag_Ltwonorm(paras_sg, list);
            % construct the detail index set ---- reduced margin
            X_diff = stochcol_rmargin_set(X, stochcol_margin_set(X));
%            X_diff = stochcol_margin_set(X);
            % Attributes of the general sparse grid interpolation corresponding
            % to the new index set, which is the union of the original index
            % set and the detail index set
            paras_sg_full = stochcol_sg([X; X_diff], rule_id);
            % set grid due to the inclusion of the detail index set
            % the grid is used to label the collocation nodes
            [gridd_diff, grid_diff_ind] = setdiff(paras_sg_full{4}, ...
                paras_sg{4}, 'rows');
            coords_diff = paras_sg_full{9}(grid_diff_ind,:);
        else
            fprintf('Spatial refinement step...\n');
            % Mesh refinement
            fprintf('original number of elements is %g\n',length(paras_fem{2}(:,1)));
            paras_fem = stochcol_mesh_refine(MMele, MMedge, paras_fem, ...
                                             paras_detail, pmethod);
            fprintf('     new number of elements is %g\n\n',length(paras_fem{2}(:,1)));
            % Collection of edge information for flux jump computation and
            % outputs of linear detail space Y grid-generator based on
            % refined mesh
            paras_detail = stochcol_mesh_detail(paras_fem);
            sols = zeros(length(paras_fem{1}), size(coords, 1));
            errests = zeros(1,size(coords, 1));
            MMeles  = cell(1,size(coords, 1));
            MMedges  = cell(1,size(coords, 1));
            for k = 1:size(coords, 1)
                [x_gal, ~] = stochcol_fem_solver(coords(k, :), ...
                    paras_fem, aa, rhs_fun);
                [errest, ~, MMele, MMedge] = stochcol_fem_estimator(pmethod, ...
                    paras_fem_errest, paras_fem, paras_detail, x_gal, coords(k, :), ...
                    aa, aax1, aax2, rhs_fun);
                sols(:, k) = x_gal;
                errests(k) = errest;
                MMeles{k} = MMele;
                MMedges{k} = MMedge;
            end
        end
    end %------------- of adaptive loop
    

    % spatial error indicator (assembled from components)
    serrest = errests*L_two_norm;
    % generate union of marked elements and edges (single grid code)
    MMele  = munion(MMeles);
    MMedge = munion(MMedges);

    % parametric error indicator
    % compute interpolated FE approximations corresponding to the new
    % collocation nodes
    sols_diff = zeros(length(paras_fem{1}), size(coords_diff, 1));
    errest2s = zeros(1, size(coords_diff, 1));
    a_unit = @(x1, x2) ones(size(x1));
    [A_unit,~] = stochcol_fem_setup(paras_fem{1}, paras_fem{2}, a_unit, rhs_fun);
    for k = 1:size(coords_diff, 1)
        % FE solution for a propective collocation node
        [x_gal, Anbc] = stochcol_fem_solver(coords_diff(k, :), ...
            paras_fem, aa, rhs_fun);
        sols_diff(:, k) = x_gal;
        % vector of multivariate Lagrange polynomials for the collocation node
        LL_diff = stochcol_getinterpolant_2(gridd, ...
            clincombiset, indlincombiset, coords_diff(k,:), polys);
        % interpolated FE solution at the collocation node
        uyy = sols*LL_diff;
        % H1 seminorm of the interpolation error at the collocation node
        errest2s(k) = sqrt((x_gal - uyy)' * A_unit * (x_gal - uyy));
    end
    
    % indexwise parametric error estimates based on the detail collocation nodes
    perrests = stochcol_est_parametric(X_diff, errest2s, ...
        gridd_diff, list, rule_id);
    % parametric error indicator
    perrest = sum(perrests);
    % marked index set from the detail index set
    [Mset, ~] = dorfler_marking(X_diff, perrests, pmthreshold);
 
    % total error indicator
    tot_err_est = serrest + perrest;
                           
    % compute direct estimates of error
    % spatial
    serrest_d = stochcol_direct_estimatorX(X, ...
        paras_sg, paras_fem, list, pmethod, rhs_fun, aa);
    % parametric
    paras_sg_diff = stochcol_sg(X_diff, rule_id);
    G_diff = stochcol_gmatrices(paras_sg_diff{4}, ...
            paras_sg_diff{5},paras_sg_diff{6}, list);
    ncpts = length(G_diff(:,1));
    ncptsf = size(coords, 1)+ size(coords_diff, 1);
    grid_old_ind = setdiff([1:ncptsf]',grid_diff_ind);
    % construct array containing all the cpoint solution vectors
    sols_all = nan(length(paras_fem{1}), ncptsf);
    sols_all(:,grid_old_ind) = sols; sols_all(:,grid_diff_ind) = sols_diff;
        if ncptsf > ncpts
        fprintf('\nredundant sparse grid solutions ..\n')
        [gdiff,indx]=setdiff(paras_sg_full{4}, paras_sg_diff{4}, 'rows');
%        gdiff
        iactive=setdiff([1:ncptsf]',indx);
        sols_all=sols_all(:,iactive);
        elseif ncptsf < ncpts, error('Oops ..  fatal logic issue'),
        end
    % compute the energy norm of the SG solution via hierarchical surplus
    perrest_d = sqrt(sum(dot(G_diff, sols_all' * A_unit * sols_all)));
    tot_err_direct = serrest_d + perrest_d;
                             
    % output error estimates
     if iter==0, fprintf('\n\nIteration %i \n',iter), end
    fprintf(['   spatial error estimate', ...
    ' is %10.4e  vs  %10.4e (spatial indicator) \n'],serrest_d,serrest)
    fprintf(['parametric error estimate', ...
    ' is %10.4e  vs  %10.4e (parametric indicator)\n'],perrest_d,perrest)
    fprintf(['overall estimate from indicators', ...
        ' is %10.4e '],tot_err_est)
    fprintf(['\n   overall direct error estimate', ...
    ' is <strong> %10.4e </strong>\n'],tot_err_direct)
    
    iter = iter + 1;
    
    err_p_iter(iter) = perrest;
    err_s_iter(iter) = serrest;
    error_iter(iter) = tot_err_est;
    
    err_p_d_iter(iter) = perrest_d;
    err_s_d_iter(iter) = serrest_d;
    error_d_iter(iter) = tot_err_direct;

    sols_iter{iter} = sols;
    paras_sg_iter{iter} = paras_sg;
    paras_fem_iter{iter} = paras_fem;
    dof(iter) = size(sols,1)*size(sols,2);

end

endLoopTime = toc(startLoopTime);

% plot error estimates
figure(98);
loglog(dof, error_iter, 'o-k', dof, error_d_iter, 's-b')
hold on, grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('$\bar\eta$ (total)', '$\eta$ (total)', ...
       'Location', 'Best','interpreter','latex');
axis tight
figure(99);
loglog(dof, error_d_iter, 'o-k', dof, err_s_d_iter, 's-b', dof, err_p_d_iter, 'd-r')
hold on
grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('\eta (total)', '\mu (spatial)', '\tau (parametric)', ...
       'Location', 'Best')
axis tight
%
figure(100);
loglog(dof, error_iter, 'o-k', dof, err_s_iter, 's-b', dof, err_p_iter, 'd-r')
hold on
grid on
xlabel('degrees of freedom')
ylabel('estimated error')
legend('$\bar\eta$ (total)', '$\bar\mu$ (spatial)', '$\bar\tau$ (parametric)', ...
       'Location', 'Best','interpreter','latex');
axis tight

% final mesh
if pmethod == 1 % P1
    plot_mesh(paras_fem{2}, paras_fem{1}, 'Final FE mesh');
    pause(1)
elseif pmethod == 2
    plot_mesh(paras_fem{6}, paras_fem{5}, 'Final FE mesh');
end
% final nodes
coords = paras_sg{9};
figure(988)
scatter(coords(:,1),coords(:,2),100,'o','filled');
axis square
box on
grid on
title('final collocation nodes in first two directions')

%% postprocess to generate statistics
multiweights = paras_sg{8};
MEAN = 0;
VAR  = 0;
for k = 1:size(coords,1)
    % total pdf
    pyy = 1;
    for m = 1:M
        pyy = pyy * fun_p(coords(k,m));
    end
    MEAN = MEAN + pyy * multiweights(k) * sols(:,k); % mean
    VAR  = VAR  + pyy * multiweights(k) * sols(:,k).^2; % variance
end
VAR = VAR - MEAN.^2;
% VAR could be negative if number of parametric pts is not large enough
if min(VAR) < 0
    error('number of collocation points is not large enough!')
end
STD = VAR.^0.5;

fprintf('\n Final sparse grid\n')
disp([paras_sg{4}])

fprintf('\n Tolerance reached in %g iterations',iter-1)
fprintf('\n    after %g parametric refinements',glevel)
fprintf('\n              Mean maximum %9.6f\n',max(MEAN))
fprintf('          Variance maximum %9.6f\n',max(VAR))
fprintf('Standard Deviation maximum %9.6f\n',max(STD))
fprintf('\n<strong>Total elapsed time:</strong> %.2f sec\n',endLoopTime);

if pmethod == 1
% plot solutions
sc_plotsol_p1(dom_type,MEAN,VAR,paras_fem{2},paras_fem{1})
end
