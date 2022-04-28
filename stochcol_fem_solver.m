function [x_gal, Anbc] = stochcol_fem_solver(coord, paras_fem, a_fun, rhs)
%STOCHCOL_FEM_SOLVER FEM solver for sample diffusion problem
%
%  [x_gal, Anbc] = stochcol_fem_solver(coord, paras_fem, a_fun, rhs)
%
%   TIFISS function: FX 28 November 2019; DJS June 2020
% Copyright (c) 2019 F. Xu & D.J. Silvester
xy     = paras_fem{1};
evt    = paras_fem{2};
bound  = paras_fem{3};
nvtx = length(xy);
% solvingTime = tic;
ayy = @(x1, x2) a_fun(x1, x2, coord);
% check the minimum of a
a_val = ayy(xy(:, 1), xy(:, 2));
% fprintf('\nmininum of diffusion coefficient %9.5f\n', min(a_val))
if min(a_val) < 0
    error('Negative diffusion coefficient!')
end
% assembly finite element matrices/vectors
[Anbc,bnbc] = stochcol_fem_setup(xy, evt, ayy, rhs);
% impose boundary condition and solve
nodes = (1:nvtx)';
intern = nodes(~ismember(nodes,bound));
xbd = xy(bound,1); % x coordinate of boundary nodes
ybd = xy(bound,2); % y coordinate of boundary nodes
bc = specific_bc(xbd,ybd);
[A,b] = imposebc(Anbc,bnbc,intern,bound,bc);
x_gal = zeros(nvtx,1);
x_gal(bound) = bc;
x_gal(intern) = A\b;
%x_gal(intern) = minres(A,b,1e-8,9999);
%fprintf('solution maximum %9.5f\n',max(x_gal))
%energy = sqrt(x_gal' * Anbc * x_gal);
%fprintf(' solution energy %9.5f\n',energy);
% fprintf('Solving took %.5f sec\n',toc(solvingTime));
