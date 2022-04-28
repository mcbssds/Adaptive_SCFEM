function [edgeres] = edgeres_p1_with_p1(xy, evt, eboundt, p1sol, eex, ...
    tve, els, coeff_fun)
%EDGERES_P1_WITH_P1 edge residuals for P1 solution using mid-edge P1 functions
%
%   [edgeres] = edgeres_p1_with_p1(xy, evt, eboundt, p1sol, eex, tve, ...
%               els, diff_fun)
%
%   input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%          eboundt    element boundary mapping matrix
%            p1sol    vertex solution vector
%              eex    element connectivity array
%              tve    edge location array
%              els    elementwise edge lengths
%        coeff_fun    function handle for the diffusion coefficient
%
%   output:
%          edgeres    edge residuals
%
% Function(s) called: gausspoints_oned
%                     tderiv
%                     p1fluxjmps
%
% See also EDGERES_P1_WITH_P2
%
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

x = xy(:,1);
y = xy(:,2);
nel = size(evt,1);

% Construct the 1D integration rule
ngpt = 7;
[oneg,onew] = gausspoints_oned(ngpt);

% Initialisation
edgeres = zeros(nel,3);

xl_v = zeros(nel,3);
yl_v = zeros(nel,3);
% Recover local coordinates
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
end

%fprintf('computing P1 flux jumps... ')

% Loop over Gaussian points
for igpt = 1:ngpt
    
    sigpt = oneg(igpt);
    sigpt_ref_edge = (1.0 + sigpt)/2.0; % [-1,1] -> [0,1]    reference-edge map
    sigpt_l = (1.0 + sigpt)/4.0;        % [-1,1] -> [0,1/2]   LEFT sub-edge map
    sigpt_r = (3.0 + sigpt)/4.0;        % [-1,1] -> [1/2,1]  RIGHT sub-edge map
    wigpt = onew(igpt);
    
    % First edge
    [~,~,phi_v_1,~,~] = tderiv(sigpt_ref_edge,1-sigpt_ref_edge,xl_v,yl_v);
    
    % Second edge
    [~,~,phi_v_2,~,~] = tderiv(0,sigpt_ref_edge,xl_v,yl_v);
    
    % Third edge
    [~,~,phi_v_3,~,~] = tderiv(sigpt_ref_edge,0,xl_v,yl_v);
    
    % Jump of the finite element solution over the left sub-element edges
    [njmp_l] = p1fluxjmps(p1sol, eex, xy, evt, eboundt, tve, sigpt_l, ...
        coeff_fun);
    
    % Jump of the finite element solution over the right sub-element edges
    [njmp_r] = p1fluxjmps(p1sol, eex, xy, evt, eboundt, tve, sigpt_r, ...
        coeff_fun);
    
    % Contribution of the first edge
    edgeres(:,1) = edgeres(:,1) + (1/2) * wigpt * njmp_l(:,1) .* phi_v_1(:,2) .* (els(:,1)/4);
    edgeres(:,1) = edgeres(:,1) + (1/2) * wigpt * njmp_r(:,1) .* phi_v_1(:,3) .* (els(:,1)/4);
    
    % Contribution of the second edge
    edgeres(:,2) = edgeres(:,2) + (1/2) * wigpt * njmp_l(:,2) .* phi_v_2(:,3) .* (els(:,2)/4);
    edgeres(:,2) = edgeres(:,2) + (1/2) * wigpt * njmp_r(:,2) .* phi_v_2(:,1) .* (els(:,2)/4);
    
    % Contribution of the third edge
    edgeres(:,3) = edgeres(:,3) + (1/2) * wigpt * njmp_l(:,3) .* phi_v_3(:,2) .* (els(:,3)/4);
    edgeres(:,3) = edgeres(:,3) + (1/2) * wigpt * njmp_r(:,3) .* phi_v_3(:,1) .* (els(:,3)/4);
    
end
% end loop over Gaussian points

%fprintf('done\n')

end  % end function
