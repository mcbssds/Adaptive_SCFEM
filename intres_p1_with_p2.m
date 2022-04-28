function [intres] = intres_p1_with_p2(xy, evt, p1sol, gradcoeffx_fun, ...
    gradcoeffy_fun, rhs_fun)
%INTRES_P1_WITH_P2 interior residuals for P1 solution using P2 bubble functions
%
%   [intres] = intres_p1_with_p2(xy,xl_s,yl_s,evt,p1sol)
%
%   input:
%               xy    vertex coordinate vector
%              evt    element mapping matrix
%            p1sol    vertex solution vector
%   gradcoeffx_fun    function handle for the first component of the
%                     gradient of the coefficient
%   gradcoeffy_fun    function handle for the second component of the
%                     gradient of the coefficient
%          rhs_fun    function handle for the right-hand side
%
%   output:
%          intres     interior residuals
%
% Function(s) called: triangular_gausspoints
%                     tgauss
% 
% See also INTRES_P1_WITH_P1
%
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  x = xy(:,1); 
  y = xy(:,2);
  nel = length(evt(:,1));

% Recover local coordinates and local solution   
  for ivtx = 1:3 
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx));
      sl_v(:,ivtx) = p1sol(evt(:,ivtx));
  end

% Construct 2D gaussian rule over the reference triangle
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);

% Preallocate matrices
  intres = zeros(nel,4);                                                      
  bde = zeros(nel,4,3);
  fde = zeros(nel,4);
 
% Loop over Gauss points 
  for igpt = 1:nngpt         
      sigpt = s(igpt);
      tigpt = t(igpt);
      wght = wt(igpt);       
          
      % Evaluate derivatives
      [jac_v,~,~,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_v,yl_v);
      [psi_v,~,~] = tqderiv(sigpt,tigpt,xl_v,yl_v);
          
      % Gradient of the diffusion coefficients
      [diffx] = tgauss(sigpt, tigpt, xl_v, yl_v, gradcoeffx_fun);
      [diffy] = tgauss(sigpt, tigpt, xl_v, yl_v, gradcoeffy_fun);

      % Source f
      [rhs_f] = tgauss(sigpt, tigpt, xl_v, yl_v, rhs_fun);
     
      % Loop over the four bubble functions
      for j = 1:4
             
          % Compute rhs-contribution from the source f
          fde(:,j) = fde(:,j) + wght * rhs_f(:) .* psi_v(:,j+3) .* jac_v(:);             
            
          % Compute div(a*grad)-contribution = grad(a)*grad(u_tau): loop
          % over vertices hat functions
          for i = 1:3
              bde(:,j,i) = bde(:,j,i) + wght * diffx(:) .* dphidx_v(:,i) .* psi_v(:,j+3);
              bde(:,j,i) = bde(:,j,i) + wght * diffy(:) .* dphidy_v(:,i) .* psi_v(:,j+3);
          end
          % end vertices hat functions loop         
      end
      % end four bubble functions loop
  end    
% end Gauss points loop
    
% Assemble interior residuals from rhs-contributions
  for i = 1:4
      intres(:,i) = intres(:,i) + fde(:,i); 
  end
  
% Multiply div(a*grad)-contribution by Galerkin solution
  for j = 1:4
      for k = 1:3
          intres(:,j) = intres(:,j) + bde(:,j,k).*sl_v(:,k);
      end
  end
  
end  % end function
