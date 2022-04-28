function [intres] = intres_p1_with_p1(xy, xl_s, yl_s, evt, p1sol, ...
    subdivPar, gradcoeffx_fun, gradcoeffy_fun, rhs_fun)
%INTRES_P1_WITH_P1 interior residuals for P1 solution using mid-edge P1 functions
%
%   [intres] = intres_p1_with_p1(xy, xl_s, yl_s, evt, p1sol, subdivPar, ...
%              gradcoeffx_fun, gradcoeffy_fun, rhs_fun)
%
%   input:
%               xy    vertex coordinate vector
%             xl_s    x-coordinates physical sub-elements
%             yl_s    y-coordinates physical sub-elements
%              evt    element mapping matrix
%            p1sol    vertex solution vector
%        subdivPar    red/bisec3 uniform sub-division switch
%   gradcoeffx_fun    function handle for the first component of the
%                     gradient of the coefficient
%   gradcoeffy_fun    function handle for the second component of the
%                     gradient of the coefficient
%          rhs_fun    function handle for the right-hand side
%   output:
%           intres    interior residual
%
% Function(s) called: triangular_gausspoints
%                     subelt_transf
%                     tderiv
%                     tgauss
%
% See also INTRES_P1_WITH_P2
%
%   TIFISS function: LR; 05 October 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi

x = xy(:,1);
y = xy(:,2);
nel = size(evt,1);

% Construct the integration rule
nngpt = 7;
[s,t,wt] = triangular_gausspoints(nngpt);

% Recover local coordinates and local solution
xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
sl_v = zeros(nel,3);
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
    sl_v(:,ivtx) = p1sol(evt(:,ivtx));
end

% Preallocate matrices
intres = zeros(nel,3);
%
bdem = zeros(nel,4,3,3);
bde = zeros(nel,3,3);
%
fdem = zeros(nel,4,3);
fde = zeros(nel,3);
%
xl_m = zeros(nel,3);
yl_m = zeros(nel,3);

% Loop over sub-elements
for subelt = 1:4
    % Recover local coordinates over sub-element
    for ivtx = 1:3
        xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
        yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
    end
    % Loop over Gauss points
    for igpt = 1:nngpt
        sigpt = s(igpt);
        tigpt = t(igpt);
        wght = wt(igpt);
        [sigptloc,tigptloc] = subelt_transf(sigpt,tigpt,subelt,subdivPar);
        
        % Evaluate derivatives, gradient of the coefficients, and rhs
        [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigptloc,tigptloc,xl_v,yl_v);
        [jac_m, ~, phi_m, ~, ~] = tderiv(sigpt, tigpt, xl_m, yl_m);
        [diffx] = tgauss(sigpt, tigpt, xl_m, yl_m, gradcoeffx_fun);
        [diffy] = tgauss(sigpt, tigpt, xl_m, yl_m, gradcoeffy_fun);
        [rhs_m] = tgauss(sigpt, tigpt, xl_m, yl_m, rhs_fun);
        
        % Loop over mid-edge hat functions
        for j = 1:3
            
            % Compute rhs-contribution from the source f
            fdem(:,subelt,j) = fdem(:,subelt,j) + wght * rhs_m(:).*phi_m(:,j).*jac_m(:);
            
            % Compute div(a*grad)-contribution = grad(a)*grad(u_h): loop
            % over vertices old hat functions
            for i = 1:3
                bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffx(:) .* dphidx_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
                bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * diffy(:) .* dphidy_v(:,i) .* phi_m(:,j) .* invjac_v(:) .* jac_m(:);
            end
            % end vertices old hat functions loop
        end
        % end mid-edge hat functions loop
    end
    % end Gauss points loop
end
% end sub-elements loop

% -----------------------------------------------------------------------------
% Manual assembly of subelement contributions
% -----------------------------------------------------------------------------
if subdivPar == 1
    %
    % Red sub-division: assembling
    %
    % First edge
    bde(:,1,1) = bdem(:,2,3,1) + bdem(:,3,2,1) + bdem(:,4,1,1);
    bde(:,1,2) = bdem(:,2,3,2) + bdem(:,3,2,2) + bdem(:,4,1,2);
    bde(:,1,3) = bdem(:,2,3,3) + bdem(:,3,2,3) + bdem(:,4,1,3);
    fde(:,1)   = fdem(:,2,3)   + fdem(:,3,2)   + fdem(:,4,1);
    % Second edge
    bde(:,2,1) = bdem(:,1,3,1) + bdem(:,3,1,1) + bdem(:,4,2,1);
    bde(:,2,2) = bdem(:,1,3,2) + bdem(:,3,1,2) + bdem(:,4,2,2);
    bde(:,2,3) = bdem(:,1,3,3) + bdem(:,3,1,3) + bdem(:,4,2,3);
    fde(:,2)   = fdem(:,1,3)   + fdem(:,3,1)   + fdem(:,4,2);
    % Third edge
    bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,1,1) + bdem(:,4,3,1);
    bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,1,2) + bdem(:,4,3,2);
    bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,1,3) + bdem(:,4,3,3);
    fde(:,3)   = fdem(:,1,2)   + fdem(:,2,1)   + fdem(:,4,3);
    
else
    %
    % Bisec3 sub-division: assembling
    %
    % First edge
    bde(:,1,1) = bdem(:,3,2,1) + bdem(:,4,2,1);
    bde(:,1,2) = bdem(:,3,2,2) + bdem(:,4,2,2);
    bde(:,1,3) = bdem(:,3,2,3) + bdem(:,4,2,3);
    fde(:,1)   = fdem(:,3,2)   + fdem(:,4,2);
    % Second edge
    bde(:,2,1) = bdem(:,1,3,1) + bdem(:,2,1,1) + bdem(:,3,3,1) + bdem(:,4,1,1);
    bde(:,2,2) = bdem(:,1,3,2) + bdem(:,2,1,2) + bdem(:,3,3,2) + bdem(:,4,1,2);
    bde(:,2,3) = bdem(:,1,3,3) + bdem(:,2,1,3) + bdem(:,3,3,3) + bdem(:,4,1,3);
    fde(:,2)   = fdem(:,1,3)   + fdem(:,2,1)   + fdem(:,3,3)   + fdem(:,4,1);
    % Third edge
    bde(:,3,1) = bdem(:,1,2,1) + bdem(:,2,2,1);
    bde(:,3,2) = bdem(:,1,2,2) + bdem(:,2,2,2);
    bde(:,3,3) = bdem(:,1,2,3) + bdem(:,2,2,3);
    fde(:,3)   = fdem(:,1,2)   + fdem(:,2,2);
end

% Assemble interior residuals from rhs source contribution
for i = 1:3
    intres(:,i) = intres(:,i) + fde(:,i);
end

% Assemble by adding the div(a*grad)-contribution times Galerkin solution
for j = 1:3
    for k = 1:3
        intres(:,j) = intres(:,j) + bde(:,j,k).*sl_v(:,k);
    end
end

end  % end function