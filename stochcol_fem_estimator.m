function [errest, elerr, MMele, MMedge] = stochcol_fem_estimator(...
    pmethod, paras_fem_errest, paras_fem, paras_detail, x_gal, coord, ...
    aa, aax1, aax2, rhs_fun, tout)
%STOCHCOL_FEM_ESTIMATOR energy error estimator
%
%  [errest, elerr, MMele, MMedge] = stochcol_fem_estimator(pmethod, ...
%    paras_fem_errest, paras_fem, paras_detail, x_gal, coord, ...
%    aa, aax1, aax2, rhs_fun, tout)
%
%   TIFISS function: FX 28 November 2019; DJS 31 January 2021
% Copyright (c) 2019 F. Xu

if nargin < 11, tout=0; end
xy      = paras_fem{1};
evt     = paras_fem{2};
eboundt = paras_fem{4};
eex     = paras_detail{1};
tve     = paras_detail{2};
els     = paras_detail{3};
evtY    = paras_detail{4};
xyY     = paras_detail{5};
boundY  = paras_detail{6};
Ybasis  = paras_detail{7};
coeff_fun = @(x1, x2) aa(x1, x2, coord);
gradcoeffx1_fun = @(x1, x2) aax1(x1, x2, coord);
gradcoeffx2_fun = @(x1, x2) aax2(x1, x2, coord);
subdivPar = paras_fem_errest(1);
% ESTIMATE
% -------------------------------------------------------------------
if tout,
fprintf('\nA posteriori error estimation\n'); end
%
% Compute a posteriori error estimation
errorTime = tic;
%
if pmethod == 1
    % P1-error estimation
    pestim = paras_fem_errest(2);
    if pestim == 1
        estimtype = paras_fem_errest(3);
        markedgelem = paras_fem_errest(4);
        markstrat = paras_fem_errest(5);
        threshold = paras_fem_errest(6);
        % Using linear midpoint Hat functions
        if tout,
        fprintf('Error estimation using 3 edge midpoint linear functions\n'); end
        if estimtype == 1
            % Hierarchical eY estimator: elementwise residual problem
            if tout,
            fprintf('Hierarchical eY estimator: solving elementwise residual problems\n'); end
            [elerr, fe, ae] = diffpost_p1_with_p1(xy, evt, eex, tve, ...
                els, eboundt, x_gal, subdivPar, coeff_fun, ...
                gradcoeffx1_fun, gradcoeffx2_fun, rhs_fun);
            [~, elerr] = diffpost_p1_bc(ae, fe, elerr, xy, evt, eboundt);
            %
            % Global error estimate
            errest = norm(elerr,2);
        elseif estimtype == 2
            % Hierarchical eY estimator: solving the assembled linear system
            if tout
            fprintf('Hierarchical eY estimator: solving assembled linear system\n'); end
            [elerr, ederr, errest] = diffpost_p1_with_p1_linsys(evt, xy, ...
                eboundt, x_gal, evtY, xyY, boundY, Ybasis, subdivPar, ...
                coeff_fun, rhs_fun);
        else% estimtype == 3
            % 2-level error estimator
            if tout
            fprintf('Two-level estimator\n'); end
            [elerr, ederr, errest] = diffpost_p1_with_p1_2level(evt, xy, ...
                eboundt, x_gal, evtY, xyY, boundY, Ybasis, subdivPar, ...
                coeff_fun, rhs_fun);
        end
    elseif pestim == 2
        markedgelem = paras_fem_errest(3);
        markstrat = paras_fem_errest(4);
        threshold = paras_fem_errest(5);
        % Using quadratic Midpoint Bubble functions
        if tout,
        fprintf('Error estimation using 4 quadratic bubble functions\n'); end
        [elerr, fe, ae] = diffpost_p1_with_p2(xy, evt, eex, tve, els, ...
            eboundt, x_gal, coeff_fun, gradcoeffx1_fun, ...
            gradcoeffx2_fun, rhs_fun);
        [~, elerr] = diffpost_p1_bc(ae, fe, elerr, xy, evt, eboundt);
        %
        % Global error estimate
        errest = norm(elerr,2);
    end
else% pmethod == 2
    markedgelem = paras_fem_errest(2);
    markstrat = paras_fem_errest(3);
    threshold = paras_fem_errest(4);
    % P2-error estimation
    if tout,
    fprintf('Error estimation using quartic bubble functions\n'); end
    %
    % check performance with that of legacy code
    %tic; [elerr_x_p4,error_total_p4,fe_p4,ae_p4] = diffpost_p2_with_p4_x(xy,evt,eboundt,x_gal); toc
    [elerr, ~, ~] = diffpost_p2_with_p4(xy, evt, eex, tve, els, ...
        eboundt, x_gal, coeff_fun, rhs_fun);
    %
    % Global error estimate
    errest = norm(elerr,2);
end % end if
%
if tout,
fprintf('Estimated energy error: %10.4e\n',errest);
fprintf('Estimation took %.5f sec\n',toc(errorTime)); end
% -------------------------------------------------------------------
% MARK
% -------------------------------------------------------------------
if markedgelem == 1
    % Marking elements
    [Mset] = marking_strategy_fa(elerr, markstrat, threshold);
else%markedgelem == 2
    % Marking edges
    [Mset] = marking_strategy_fa(ederr, markstrat, threshold);
end
%
% Overall set of marked elements (and edges)
[MMele, MMedge] = get_all_marked_elem(Ybasis, evtY, Mset, markedgelem);
