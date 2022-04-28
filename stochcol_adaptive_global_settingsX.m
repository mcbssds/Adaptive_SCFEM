%----------------------------------------- legacy code
clear
close ALL
clc
%% set parameters
fprintf('choose type of diffusion coefficient ');
fprintf('\n     1.  exponential of truncated affine expansion')
fprintf('\n     2.  square of truncated affine expansion');
fprintf('\n     3.  standard truncated affine (K-L) expansion\n')
rf_type = default('default is exponential case',1);
if ~ismember(rf_type,[1, 2, 3])
    error('Wrong selection for type of diffusion coefficient!')
end

fprintf('\nchoose type of spatial expansion coefficient ');
fprintf('\n     1.  separable exponential ')
fprintf('\n     2.  synthetic (Eigel slow) expansion\n');
expansion_type  = default('default is separable exponential',1);
if ~ismember(expansion_type,[1,2])
    error('Wrong selection for expansion type!')
end
if expansion_type == 1 % seperable exponential
    sigma = default('SE standard deviation (default is 0.5)', 0.5);
    ax = 1;
    ay = 1;
    correl_x = default('correlation length in x (default is 1)', 1);
    correl_y = default('correlation length in y (default is 1)', 1);
elseif expansion_type == 2 % Eigel
    sigma = default('Eigel standard deviation (default is 0.547)', 0.547);
end

fprintf('\nchoose type of random variable ');
fprintf('\n     1.  uniform ')
fprintf('\n     2.  truncated Gaussian\n');
rv_id = default('default is uniform', 1);
if ~ismember(rv_id,[1,2])
    error('Wrong selection for type of random variable!')
end
if rv_id == 2
    sigma = 0.5;
end

fprintf('\nchoose type of finite element approximation');
fprintf('\n     1.  P1 ')
fprintf('\n     2.  P2\n');
pmethod = default('default is P1', 1);
if ~ismember(pmethod,[1,2])
    error('Wrong selection for approximation method!')
end

% Red/Bisec3 for spatial error estimation 1/2? See:
% Ainsworth, Oden, A posteriori error estimation in finite element analysis,
% Wiley, 2000 - Figure 5.2 (p. 87) for the basis functions in both cases.
subdivPar = 2;

% Error estimation and estimators type
paras_fem_errest = subdivPar;
if pmethod == 1
    pestim = default('\nError estimation: linear/quadratic bubble functions 1/2? (default 1)',1);
    paras_fem_errest = [paras_fem_errest, pestim];
    if pestim == 1
        % Error estimation type (for P1 approximations only)
        % 1 - eY hierarchical estimator (elementwise residual problems)
        % 2 - eY hierarchical estimator (assembled system for the residual problems)
        % 3 - 2-level error estimator
        fprintf('Estimator type:\n');
        fprintf('   1. hierarchical estimator (elementwise residual problems)\n');
        fprintf('   2. hierarchical estimator (fully assembled system for residual problem)\n');
        fprintf('   3. 2-level estimator\n');
        estimtype = default('(default 1)',1);
        paras_fem_errest = [paras_fem_errest, estimtype];
        if ~ismember(estimtype,[1,2,3]), error('Estimator type not allowed!'); end
        %
        % Marking elements or edges 1/2?
        % This depends on the estimator type:
        % - ypestim=1   -> only elements can be marked, i.e., markedgelem=1;
        % - ypestim=2/3 -> both elements and edges can be marked, i.e.,markedgelem=1/2.
        if estimtype == 1
            % Marking elements only
            markedgelem = 1;
        elseif estimtype == 2
            markedgelem = default('Marking elements/edges 1/2 (default 1)',1);
        else%estimtype = 3
            markedgelem = default('Marking elements/edges 1/2 (default 2)',2);
        end
        paras_fem_errest = [paras_fem_errest, markedgelem];
        if ~ismember(markedgelem,[1,2])
            error('Marking type not allowed!');
        end
    elseif pestim == 2
        % Marking elements
        fprintf('Using hierarchical estimator (elementwise residual problems)\n');
        markedgelem = 1;
        paras_fem_errest = [paras_fem_errest, markedgelem];
    else
        error('Estimation type not allowed!');
    end
elseif pmethod == 2
    % Marking elements
    fprintf('Using hierarchical estimator (elementwise residual problems)\n');
    markedgelem = 1;
    paras_fem_errest = [paras_fem_errest, markedgelem];
end

% Marking threshold parameters for both elements/edges and indices:
% 1 - maximum strategy:
%     large threshold -> small set of marked elements/edges
%     small threshold -> large set of marked elements/edges
% 2 - Doerfler (equilibration) strategy:
%     large threshold -> large set of marked elements/edges
%     small threshold -> small set of marked elements/edges
markstrat   = default('Marking strategy: maximum or equilibration 1/2? (default 2)',2);
smthreshold = default('Threshold parameter (default 0.3)',0.3);
paras_fem_errest = [paras_fem_errest, markstrat, smthreshold];

% domain type
dom_type = default('domain type: Square/L-shaped 1/2 (default is 1)',1);
if dom_type == 1
    dom_paras = dom_type;
elseif dom_type == 2
    dom_paras = dom_type;
    mesh_type = default('\nStructured/unstructured mesh 1/2 (default 1)',1);
    dom_paras = [dom_paras, mesh_type];
end

% parameter space [-L, L]^M
L = 1;
M = default('dimension of parametric space (default is 4)',4);

fprintf('\nchoose type of collocation nodes');
fprintf('\n     1.  Leja ')
fprintf('\n     2.  CC\n');
rule_id = default('default is CC nodes',2);

% KL expansion
if expansion_type == 1 % seperable exponential
    input = [M, ax, ay, correl_x, correl_y, sigma];
elseif expansion_type == 2 % Eigel
    input = [M, sigma, 2];
end
KL_DATA = stoch_kl(input,expansion_type);
% function handle for diffusion coefficient w.r.t. spatial and parametric
% variabels
[aa, a] = stochcol_diffusion_coeff_fun(KL_DATA, rf_type, M);
% function handle for gradients of diffusion coefficient
[aax1, aax2] = stochcol_diffusion_gradcoeff_fun(KL_DATA, rf_type, M, aa, a);
% function handle for the right-hand side w.r.t. spatial variables
rhs_fun = @(x1,x2) ones(size(x1));

% probability density function
if rv_id == 1 % uniform
    fun_p = @(x) 0.5;
elseif rv_id == 2 % truncated Gaussian
    sigma = 1;
    fun_p = @(x) ...
        exp(-x.^2./sigma./sigma./2)...
        /sigma/sqrt(2*pi)/erf(L/sqrt(2)/sigma);
end

% 1D Lagrange polynomials
if rule_id == 1
    max_level = 9;
elseif rule_id == 2
    max_level = 7;
end
if exist('precomputation.mat', 'file') == 2 % check existence of the file precomputation.mat
    load precomputation.mat
    if max_level_p == max_level && rule_id_p == rule_id && rv_id_p == rv_id && L_p == L
    else
        polys = stochcol_onedlagpolys(max_level, rule_id);
        list = stochcol_uni_int(fun_p, polys, L);
        max_level_p = max_level;
        rule_id_p = rule_id;
        rv_id_p = rv_id;
        L_p = L;
        save precomputation.mat max_level_p rule_id_p rv_id_p fun_p polys L_p list
    end
else
    polys = stochcol_onedlagpolys(max_level, rule_id);
    list = stochcol_uni_int(fun_p, polys, L);
    max_level_p = max_level;
    rule_id_p = rule_id;
    rv_id_p = rv_id;
    L_p = L;
    save precomputation.mat max_level_p rule_id_p rv_id_p fun_p polys L_p list
end

% Marking threshold parameters for indices
pmthreshold = default('Threshold parameter for marking indices (default 0.3)',0.3);
