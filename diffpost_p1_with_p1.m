function [elerr_p, fe, Ae] = diffpost_p1_with_p1(xy, evt, eex, tve, els, ...
    eboundt, p1sol, subdivPar, coeff_fun, gradcoeffx_fun, ...
    gradcoeffy_fun, rhs_fun)
%DIFFPOST_P1_WITH_P1 a posteriori estimation for P1 using mid-edge P1 functions
%
%   [elerr_p, fe, Ae] = diffpost_p1_with_p1(xy, evt, eex, tve, els, ...
%                       eboundt, p1sol, subdivPar, coeff_fun, ...
%                       gradcoeffx_fun, gradcoeffy_fun, rhs_fun)
%
%   input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%          eboundt    element boundary mapping matrix
%            p1sol    P1 solution for diffusion problem
%        subdivPar    red/bisec3 uniform sub-division switch
%        coeff_fun    function handle for the diffusion coefficient
%   gradcoeffx_fun    function handle for the first component of the
%                     gradient of the coefficient
%   gradcoeffy_fun    function handle for the second component of the
%                     gradient of the coefficient
%          rhs_fun    function handle for the right-hand side
%
%   output:
%          elerr_p    elementwise error estimate
%               fe    elementwise rhs vectors
%               Ae    elementwise Poisson problem matrices
%
% Function(s) called: triangular_gausspoints
%                     tgauss
%                     intres_p1_with_p1
%                     edgeres_p1_with_p1
%
% See also DIFFPOST_P1_WITH_P2
%
%   TIFISS function: LR; 05 October 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi

% NOTE:
% -------------------------------------------------------------------------
% The function allows the error estimation based on mesh subdivision in
% sub-elements for both uniform (red) and criss-cross subdivisions. They
% are both based on the same edge midpoints but the corresponding linear
% basis functions will have a different support.
%
% > The parameter 'subdivPar' decide the type of subdivision
%
% You can have a look at fig. 5.2 (pag 87) of the book:
%   Ainsworth, Oden - 'A posteriori error estimation in finite element
%   analysys', Wiley, New York, 2000
% to see how the basis functions look like on uniform and criss-cross mesh

%fprintf('P1 local error estimator using mid-edge functions...\n');
x = xy(:,1);
y = xy(:,2);
nel = size(evt,1);
elerr_p = zeros(nel,1);

% Recover local coordinates
xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
end

% Construct the integration rule
nngpt = 7;
[s,t,wt] = triangular_gausspoints(nngpt);

% LHS of the linear system
% -------------------------------------------------------------------------
nnode = 3;

% Preallocate local matrices
adem = zeros(nel,4,nnode,nnode);
ae = zeros(nel,nnode,nnode);
%
xl_s = zeros(nel,4,3);
yl_s = zeros(nel,4,3);
%
xl_m = zeros(nel,3);
yl_m = zeros(nel,3);

% -----------------------------------------------------------------------------
% STEP 1: coordinates of midpoints: four sub-elements
% -----------------------------------------------------------------------------

% First physical mid-edge point
xedge1(:,1)= 0.5*(xl_v(:,2) + xl_v(:,3));
yedge1(:,1)= 0.5*(yl_v(:,2) + yl_v(:,3));

% Second physical mid-edge point
xedge2(:,1)= 0.5*(xl_v(:,1) + xl_v(:,3));
yedge2(:,1)= 0.5*(yl_v(:,1) + yl_v(:,3));

% Third physical mid-edge point
xedge3(:,1)= 0.5*(xl_v(:,1) + xl_v(:,2));
yedge3(:,1)= 0.5*(yl_v(:,1) + yl_v(:,2));

% Define the local sub-division
if subdivPar == 1
    %
    % Red sub-division
    %
    % First physical sub-element
    xl_s(:,1,1) = xl_v(:,1);      yl_s(:,1,1) = yl_v(:,1);
    xl_s(:,1,2) = xedge3(:);      yl_s(:,1,2) = yedge3(:);
    xl_s(:,1,3) = xedge2(:);      yl_s(:,1,3) = yedge2(:);
    % Second physical sub-element
    xl_s(:,2,1) = xedge3(:);      yl_s(:,2,1) = yedge3(:);
    xl_s(:,2,2) = xl_v(:,2);      yl_s(:,2,2) = yl_v(:,2);
    xl_s(:,2,3) = xedge1(:);      yl_s(:,2,3) = yedge1(:);
    % Third physical sub-element
    xl_s(:,3,1) = xedge2(:);      yl_s(:,3,1) = yedge2(:);
    xl_s(:,3,2) = xedge1(:);      yl_s(:,3,2) = yedge1(:);
    xl_s(:,3,3) = xl_v(:,3);      yl_s(:,3,3) = yl_v(:,3);
    % Fourth physical sub-element
    xl_s(:,4,1) = xedge1(:);      yl_s(:,4,1) = yedge1(:);
    xl_s(:,4,2) = xedge2(:);      yl_s(:,4,2) = yedge2(:);
    xl_s(:,4,3) = xedge3(:);      yl_s(:,4,3) = yedge3(:);
else
    %
    % Bisec3 sub-division
    %
    % First physical sub-element
    xl_s(:,1,1) = xl_v(:,1);      yl_s(:,1,1) = yl_v(:,1);
    xl_s(:,1,2) = xedge3(:);      yl_s(:,1,2) = yedge3(:);
    xl_s(:,1,3) = xedge2(:);      yl_s(:,1,3) = yedge2(:);
    % Second physical sub-element
    xl_s(:,2,1) = xedge2(:);      yl_s(:,2,1) = yedge2(:);
    xl_s(:,2,2) = xedge3(:);      yl_s(:,2,2) = yedge3(:);
    xl_s(:,2,3) = xl_v(:,2);      yl_s(:,2,3) = yl_v(:,2);
    % Third physical sub-element
    xl_s(:,3,1) = xl_v(:,2);      yl_s(:,3,1) = yl_v(:,2);
    xl_s(:,3,2) = xedge1(:);      yl_s(:,3,2) = yedge1(:);
    xl_s(:,3,3) = xedge2(:);      yl_s(:,3,3) = yedge2(:);
    % Fourth physical sub-element
    xl_s(:,4,1) = xedge2(:);      yl_s(:,4,1) = yedge2(:);
    xl_s(:,4,2) = xedge1(:);      yl_s(:,4,2) = yedge1(:);
    xl_s(:,4,3) = xl_v(:,3);      yl_s(:,4,3) = yl_v(:,3);
end

% -----------------------------------------------------------------------------
% STEP 2: left-Hand side of the linear system
% -----------------------------------------------------------------------------
for subelt = 1:4
    % Recover local coordinates of sub-elements
    for ivtx = 1:3
        xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
        yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
    end
    % Loop over Gauss points
    for igpt = 1:nngpt
        sigpt = s(igpt);
        tigpt = t(igpt);
        wght = wt(igpt);
        % Evaluate derivatives and coefficients
        [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigpt, tigpt, xl_m, yl_m);
        [coeff] = tgauss(sigpt, tigpt, xl_m, yl_m, coeff_fun);
        % Loop over the three mid-edge linear functions
        for j = 1:3
            for i = 1:3
                adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * coeff(:) .* dphidx_v(:,i) .* dphidx_v(:,j) .* invjac_v(:);
                adem(:,subelt,i,j) = adem(:,subelt,i,j) + wght * coeff(:) .* dphidy_v(:,i) .* dphidy_v(:,j) .* invjac_v(:);
            end
        end
        % end mid-edge linear functions loop
    end
    % end of Gauss point loop
end
% end subdivided element loop

% -----------------------------------------------------------------------------
% Manual assembly of subelement contributions
% -----------------------------------------------------------------------------
if subdivPar == 1
    %
    % Red sub-division: assembling
    %
    % First edge
    ae(:,1,1) = adem(:,2,3,3) + adem(:,3,2,2) + adem(:,4,1,1);
    ae(:,1,2) = adem(:,3,2,1) + adem(:,4,1,2);
    ae(:,1,3) = adem(:,2,3,1) + adem(:,4,1,3);
    % Second edge
    ae(:,2,1) = adem(:,3,1,2) + adem(:,4,2,1);
    ae(:,2,2) = adem(:,1,3,3) + adem(:,3,1,1) + adem(:,4,2,2);
    ae(:,2,3) = adem(:,1,3,2) + adem(:,4,2,3);
    % Third edge
    ae(:,3,1) = adem(:,2,1,3) + adem(:,4,3,1);
    ae(:,3,2) = adem(:,1,2,3) + adem(:,4,3,2);
    ae(:,3,3) = adem(:,1,2,2) + adem(:,2,1,1) + adem(:,4,3,3);
else
    %
    % Bisec3 sub-division: assembling
    %
    % First edge
    ae(:,1,1) = adem(:,3,2,2) + adem(:,4,2,2);
    ae(:,1,2) = adem(:,3,2,3) + adem(:,4,2,1);
    % ae(:,1,3) = empty
    % Second edge
    ae(:,2,1) = adem(:,3,3,2) + adem(:,4,1,2);
    ae(:,2,2) = adem(:,1,3,3) + adem(:,2,1,1) + adem(:,3,3,3) + adem(:,4,1,1);
    ae(:,2,3) = adem(:,1,3,2) + adem(:,2,1,2);
    % Third edge
    % ae(:,3,1) = empty
    ae(:,3,2) = adem(:,1,2,3) + adem(:,2,2,1);
    ae(:,3,3) = adem(:,1,2,2) + adem(:,2,2,2);
end

% Saving the LHS matrix of the system for output before factorization
Ae = ae;

% -----------------------------------------------------------------------------
% STEP 3: right-hand side of the linear system
% -----------------------------------------------------------------------------

% Element residual
[res_int] = intres_p1_with_p1(xy, xl_s, yl_s, evt, p1sol, subdivPar, ...
    gradcoeffx_fun, gradcoeffy_fun, rhs_fun);

% Edge residual: this is independent from the subdivision chosen
[res_edge] = edgeres_p1_with_p1(xy, evt, eboundt, p1sol, eex, tve, els, ...
    coeff_fun);
%fprintf('internal_res = %7.4e | edge_res = %7.4e\n',norm(res_int),norm(res_edge));

% Final RHS of the linear system
fe = res_int - res_edge;

% LDL^t factorisation of the LHS matrix
% -----------------------------------------------------------------------------
[ae] = ldlt_fact(ae,nel);

% -----------------------------------------------------------------------------
% STEP 4: solving the system (forward-backward substitutions)
% -----------------------------------------------------------------------------
xx = element_lusolve(ae,fe);
elerr = xx';

% Elementwise error estimation
for ivtx = 1:nnode
    elerr_p(:) = elerr_p(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
elerr_p = sqrt(elerr_p);

%fprintf('Estimated energy error: %10.4e\n',norm(elerr_p,2));

end  % end function


% ----------------------------------------------------------------------
% Child function
% ----------------------------------------------------------------------
function [ade] = ldlt_fact(ade,nel)
% LDLT factorization of the matrix ade

nn = 3; % number of hat functions per element
dd = zeros(nel,nn);
rr = zeros(nel,nn);

for kk = 1:nn-1
    for pp = 1:kk-1
        rr(1:nel,pp) = dd(1:nel,pp).*ade(1:nel,kk,pp);
    end
    dd(1:nel,kk) = ade(1:nel,kk,kk);
    for pp = 1:kk-1
        dd(1:nel,kk) = dd(1:nel,kk) - ade(1:nel,kk,pp).*rr(1:nel,pp);
    end
    for ii = kk+1:nn
        for pp = 1:kk-1
            ade(1:nel,ii,kk) = ade(1:nel,ii,kk) - ade(1:nel,ii,pp).*rr(1:nel,pp);
        end
        ade(1:nel,ii,kk) = ade(1:nel,ii,kk)./dd(1:nel,kk);
    end
end

for pp = 1:nn-1
    rr(1:nel,pp) = dd(1:nel,pp).*ade(1:nel,nn,pp);
end

dd(1:nel,nn) = ade(1:nel,nn,nn);

for pp = 1:nn-1
    dd(1:nel,nn) = dd(1:nel,nn) - ade(1:nel,nn,pp).*rr(1:nel,pp);
end

% Overwrite diagonal entries
for kk = 1:nn
    ade(1:nel,kk,kk) = dd(1:nel,kk);
end

end % end child function
