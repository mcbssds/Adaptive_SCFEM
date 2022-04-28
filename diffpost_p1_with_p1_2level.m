function [elerr, ederr, errest] = diffpost_p1_with_p1_2level(evt, xy, ...
    eboundt, xgal, evtY, xyY, boundY, Ybasis, subdivPar, coeff_fun, rhs_fun)
%DIFFPOST_P1_WITH_P1_2LEVEL computes 2-level error estimator for stochastic Galerkin P1 solution
%
% [elerr,ederr,errest] = diffpost_p1_with_p1_2level(evt,xy,eboundt,...
%                        x_gal,evtY,xyY,boundY,Ybasis,subdivPar)
%
% input:
%              evt    element mapping matrix
%               xy    vertex coordinate vector  
%          eboundt    element boundary mapping matrix
%             xgal    P1 solution for diffusion problem
%             evtY    element mapping matrix for midpoints
%              xyY    vertex coordinate vector for midpoints
%           boundY    boundary midpoints vector
%           Ybasis    Y-basis element-positions matrix
%        subdivPar    red/bisec3 uniform sub-division switch
%        coeff_fun    function handle for the diffusion coefficient
%          rhs_fun    function handle for the right-hand side
%
% output:
%            elerr    vector of 2-level element indicators
%            ederr    vector of 2-level edge indicators
%           errest    global 2-level error estimate
%
% The element indicators are recovered by square summing the 3 edge 
% indicators per element.
%
% Main reference for this estimator (in the stochastic framework):
% [BPRR18] Bespalov, Praetorius, Rocchi, Ruggeri, Goal-oriented error estimation 
% and adaptivity for elliptic PDEs with parametric or uncertaint inputs, 2018. 
%                    
% Function(s) called: triangular_gausspoints
%                     tderiv
%                     tgauss
%                     subelt_transf
%                     specific_bc
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri
  
  nvtx = size(xy,1);      % Number of vertices  (i.e., number of X-basis functions) 
  nyb  = size(Ybasis,1);  % Number of midpoints (i.e., number of Y-basis functions (boudary ones included))

% Extract elements and associated edges for Y-basis functions
  elems = Ybasis(:,[1,2]);   
  edges = Ybasis(:,[3,4]);  
 
% -----------------------------------------------------------------------------                 
% STEP 1: elementwise contributions from uniform red/bisec3 refinement type
% -----------------------------------------------------------------------------      
  [~, bde, ~ , fde] = diffpost_p1_detcontrib(xy, evt, subdivPar, coeff_fun, rhs_fun);
%  [ae, bde, ~ , fde] = diffpost_p1_detcontrib(xy, evt, subdivPar, coeff_fun, rhs_fun);
  
% -----------------------------------------------------------------------------                 
% STEP 2: assembling ||\grad\phi_j||^2_L^2(D), for all \phi_j in Y             
% -----------------------------------------------------------------------------
% Indices: elements, Y-basis, Y-basis
  ind1 = {elems(:,1), edges(:,1), edges(:,1)};
  ind2 = {elems(:,2), edges(:,2), edges(:,2)};
% Computing ae for the unit coefficient
  unit_coeff = @(x1, x2) ones(size(x1));
  [ae, ~, ~ , ~] = diffpost_p1_detcontrib(xy, evt, subdivPar, unit_coeff, rhs_fun);

% Summing contributions  
  B0phi = ae( sub2ind(size(ae) , ind1{:}) ) + ae( sub2ind(size(ae) , ind2{:}) );
  
% For a Y-basis on the boundary, we have to halve the contributions, since in 
% the above line we doubled the contributions coming from the boundary elements
  B0phi(boundY) = B0phi(boundY) ./ 2;
     
% -----------------------------------------------------------------------------                       
% STEP 3: assembling F(v) for all basis \phi_j in Y
% -----------------------------------------------------------------------------
% Indices: elements, Y-basis
  ind1 = {elems(:,1), edges(:,1)};
  ind2 = {elems(:,2), edges(:,2)};
% Summing contributions
  ff = fde( sub2ind(size(fde) , ind1{:}) ) + fde( sub2ind(size(fde) , ind2{:}) );
% For a Y-basis on the boundary, we have to halve the contributions, since in 
% the above line we doubled the contributions coming from the boundary elements
  ff(boundY) = ff(boundY) ./ 2; 
  
% -----------------------------------------------------------------------------                                
% STEP 4: Assembling B(xgal,v) for all basis \phi_j in Y
% -----------------------------------------------------------------------------  
% The term B(xgal,v) will be a vector nyb-by-1     
  Km = sparse(nyb,nvtx);
  for kYb = 1:3
      indY = evtY(:,kYb);	 
      for kXb = 1:3
          indX = evt(:,kXb);
          Km = Km  + sparse(indY, indX, bde(:, kYb, kXb), nyb, nvtx);
      end
  end  
  Bx = Km * xgal;
      
% Global two-level numerator (rhs): F(v) - B(xgal,v)
  rhs = ff - Bx;
   
% -----------------------------------------------------------------------------
% STEP 5: imposing interpolated error as Dirichlet bcs on Y-boundary nodes
% -----------------------------------------------------------------------------  
  [rhs,~,nonzerobndY] = nonzerobc_yspace(evt,evtY,xy,xyY,eboundt,rhs);
  
% -----------------------------------------------------------------------------
% STEP 6: compute the total 2-level estimate  
% -----------------------------------------------------------------------------

% Interior Y-nodes
  interiorY = 1:boundY(1)-1;
  
% Vector of edge indicators (internal and boundary Y-basis functions)
  ederr_sq = zeros(nyb,1);
  ederr_sq(interiorY)   = (rhs(interiorY).^2)   ./ B0phi(interiorY);
  ederr_sq(nonzerobndY) = (rhs(nonzerobndY).^2) ./ B0phi(nonzerobndY);
% NOTE that the second contribution (due to nonzerobndY) is zero if bcs = 0 everywhere;
  
% Vector of edge-indicators
  ederr = sqrt( ederr_sq );
  
% Global 2-level estimate 
  errest = sqrt( sum( ederr.^2 ) );  
  
% -----------------------------------------------------------------------------
% STEP 7: recovering the element indicators from edge indicators
% -----------------------------------------------------------------------------  
  [elerr] = get_errelem(evtY,ederr);  
     
end % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [ade, bde, bdem, fde] = diffpost_p1_detcontrib(xy, evt, ...
    subdivPar, coeff_fun, rhs_fun)
%Output:
%       ade     elementwise lhs contribution 
%       bde     elementwise rhs contribution from Galerkin solution
%      bdem     sub-elementwise rhs contribution from Galerkin solution
%       fde     elementwise rhs contribution from the source f
%
% The function computes the elementwise contributions to the spatial error
% over either the red or bisec3 uniform refinement.
%
% The following discrete formulation is considered
%
%   B0(eY,v) = F(v) - B(ugal,v) for all v \in Y
%
% and the elementwise contributions from the deterministic matrices on the 
% LHS as well as the contributions from the source f and the Galerkin 
% solution ugal on the RHS.
%
% NOTE that the function replicates the first part of the function
% DIFFPOST_P1_WITH_P1 and it is based on INTRES_P1_WITH_P1 for 
% contributions coming from the source and the Galerkin solution 

  nel = size(evt,1);    % Number of elements
   
% Construct the integration rule (3/7/19/73 Gaussian points)
  nngpt = 7;
  [s,t,wt] = triangular_gausspoints(nngpt);  
  
% Recover local coordinates  
  xl_v = zeros(nel,3); 
  yl_v = zeros(nel,3); 
  for ivtx = 1:3
      xl_v(:,ivtx) = xy(evt(:,ivtx),1);
      yl_v(:,ivtx) = xy(evt(:,ivtx),2);
  end
   
% Initialise local matrices  
  adem = zeros(nel,4,3,3);
  ade  = zeros(nel,3,3);   
  xl_s = zeros(nel,4,3); 
  yl_s = zeros(nel,4,3);
  xl_m = zeros(nel,3);
  yl_m = zeros(nel,3);
  
% ----------------------------------------------------------------------------- 
% STEP 1: coordinates of midpoints: four sub-elements
% -----------------------------------------------------------------------------
% First physical mid-edge point
  xedge1(:,1) = 0.5 * (xl_v(:,2) + xl_v(:,3));    
  yedge1(:,1) = 0.5 * (yl_v(:,2) + yl_v(:,3));
  
% Second physical mid-edge point
  xedge2(:,1) = 0.5 * (xl_v(:,3) + xl_v(:,1));  
  yedge2(:,1) = 0.5 * (yl_v(:,3) + yl_v(:,1));
  
% Third physical mid-edge point
  xedge3(:,1) = 0.5 * (xl_v(:,1) + xl_v(:,2));  
  yedge3(:,1) = 0.5 * (yl_v(:,1) + yl_v(:,2));

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
% STEP 3: integrals in B(eY,v) for all basis v \in Y
% ----------------------------------------------------------------------------- 
% Loop over sub-elements
  for subelt = 1:4 
      % Recover local coordinates of the current subelement
      for ivtx = 1:3
          xl_m(:,ivtx) = xl_s(:,subelt,ivtx);
          yl_m(:,ivtx) = yl_s(:,subelt,ivtx);
      end
      % Loop over Gauss points
      for igpt = 1:nngpt              
          sigpt = s(igpt);
          tigpt = t(igpt);
          wght  = wt(igpt);  
          % Evaluate derivatives
          [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_m,yl_m);
          % Evaluate variable diffusion coefficients
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
% end of subdivided element loop

% Manual assembly of subelement contributions
  if subdivPar == 1
      %
      % Red sub-division: assembling contributions
      % 
      % First edge
      ade(:,1,1) = adem(:,2,3,3) + adem(:,3,2,2) + adem(:,4,1,1);
      ade(:,1,2) = adem(:,3,2,1) + adem(:,4,1,2);
      ade(:,1,3) = adem(:,2,3,1) + adem(:,4,1,3);
      % Second edge
      ade(:,2,1) = adem(:,3,1,2) + adem(:,4,2,1);
      ade(:,2,2) = adem(:,1,3,3) + adem(:,3,1,1) + adem(:,4,2,2);
      ade(:,2,3) = adem(:,1,3,2) + adem(:,4,2,3);  
      % Third edge     
      ade(:,3,1) = adem(:,2,1,3) + adem(:,4,3,1);
      ade(:,3,2) = adem(:,1,2,3) + adem(:,4,3,2);
      ade(:,3,3) = adem(:,1,2,2) + adem(:,2,1,1) + adem(:,4,3,3); 
  else
      % 
      % Bisec3 sub-division: assembling contributions
      %
      % First edge
      ade(:,1,1) = adem(:,3,2,2) + adem(:,4,2,2);
      ade(:,1,2) = adem(:,3,2,3) + adem(:,4,2,1);
      % ae(:,1,3) = empty
      % Second edge
      ade(:,2,1) = adem(:,3,3,2) + adem(:,4,1,2);
      ade(:,2,2) = adem(:,1,3,3) + adem(:,2,1,1) + adem(:,3,3,3) + adem(:,4,1,1);
      ade(:,2,3) = adem(:,1,3,2) + adem(:,2,1,2);  
      % Third edge     
      % ae(:,3,1) = empty 
      ade(:,3,2) = adem(:,1,2,3) + adem(:,2,2,1);
      ade(:,3,3) = adem(:,1,2,2) + adem(:,2,2,2);
  end   
      
% ----------------------------------------------------------------------------- 
% STEP 3: integrals in F(v) and B(xgal,v) for all basis v in Y
% ----------------------------------------------------------------------------- 
% Preallocate matrices
  bdem = zeros(nel,4,3,3);
  bde  = zeros(nel,3,3);
  fdem = zeros(nel,4,3);
  fde  = zeros(nel,3);
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
          wght  = wt(igpt);       
          [sigptloc,tigptloc] = subelt_transf(sigpt,tigpt,subelt,subdivPar);  
          %
          % Evaluate derivatives, coefficient, and rhs
          [~,invjac_v,~,dphidx_v,dphidy_v]  = tderiv(sigptloc,tigptloc,xl_v,yl_v);
          [jac_m,~,phi_m,dphidx_m,dphidy_m] = tderiv(sigpt,tigpt,xl_m,yl_m); 
          [coeff] = tgauss(sigpt, tigpt, xl_m, yl_m, coeff_fun);
          [rhs_m] = tgauss(sigpt, tigpt, xl_m, yl_m, rhs_fun);
          %
          % Loop over mid-edge hat functions
          for j = 1:3  
              % Loop over X-basis functions
              % Contributions: \int_subelt a(x) \grad(Xbasis) \cdot \grad(Ybasis) dx 
              for i = 1:3  
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * coeff(:) .* dphidx_v(:,i) .* dphidx_m(:,j) .* invjac_v(:);            
                  bdem(:,subelt,j,i) = bdem(:,subelt,j,i) + wght * coeff(:) .* dphidy_v(:,i) .* dphidy_m(:,j) .* invjac_v(:);
              end
              % Contribution from the source f
              fdem(:,subelt,j) = fdem(:,subelt,j) + wght * rhs_m(:) .* phi_m(:,j) .* jac_m(:);  
          end
          % end mid-edge hat functions loop
      end
      % end Gauss points loop
  end
% end sub-elements loop

% Manual assembly of subelement contributions
  if subdivPar == 1
      %
      % Uniform sub-division: assembling contributions
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
      % Bisec3 sub-division: assembling contributions
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
 
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [Fv,zerobndY,nonzerobndY] = nonzerobc_yspace(evt,evtY,xy,xyY,eboundt,Fv)
%Set boundary conditions on edges-midpoints (Y-space)

% -----------------------------------------------------------------------------
% STEP 1: impose boundary conditions on all Y-boundary nodes
% -----------------------------------------------------------------------------

% Extract boundary elements with boundary edges respectively = 1, 2, and 3
  beled1 = eboundt( eboundt(:,2) == 1 , 1);
  beled2 = eboundt( eboundt(:,2) == 2 , 1);
  beled3 = eboundt( eboundt(:,2) == 3 , 1);
  
% For boundary edge = 1, the X-nodes positions = (2,3) and Y-node position = 1
  bnodesX1 = evt ( beled1, [2 3]);
  bnodeY1  = evtY( beled1, 1);
  [error1,bc_nodeY1] = interpolated_error_ybc(xy,xyY,bnodesX1,bnodeY1);
 
% For boundary edge = 2, the X-nodes positions = (3,1) and Y-node position = 2
  bnodesX2 = evt ( beled2, [3 1]);
  bnodeY2  = evtY( beled2, 2);
  [error2,bc_nodeY2] = interpolated_error_ybc(xy,xyY,bnodesX2,bnodeY2);

% For boundary edge = 3, the X-nodes positions = (1,2) and Y-node position = 3
  bnodesX3 = evt ( beled3, [1 2]);
  bnodeY3  = evtY( beled3, 3);
  [error3,bc_nodeY3] = interpolated_error_ybc(xy,xyY,bnodesX3,bnodeY3);
  
% Update the vector Fv in all boundary-node positions  
  Fv(bnodeY1) = error1;
  Fv(bnodeY2) = error2;
  Fv(bnodeY3) = error3;
  
% DEBUG: the following matrices contain on each row the Y-boundary node, the
% corresponding boundary condition, and the corresponding error. They are
% divided w.r.t. the Y-boundary nodes with positions 1,2, or 3, respectively.
%   [bnodeY1, bc_nodeY1, error1]
%   [bnodeY2, bc_nodeY2, error2]
%   [bnodeY3, bc_nodeY3, error3]
% Note that boundY = { bnodeY1, bnodeY2, bnodeY3 }

% -----------------------------------------------------------------------------
% Separate Y-boundary nodes with nonzero boundary conditions
% -----------------------------------------------------------------------------
  
% Among all the Y-boundary nodes, we want to return those ones for which 
% the corresponding bcs are non-homogeneous (if any).
% However, we have to ensure that where bc == 0 the returned value is really 
% zero (and not some value ~1.3e-17) because this would not be read as zero. 
% To this end, if some values is returned as ~1.3e-17, we find it by looking
% at those values smaller than, for instance, 1e-10.

% Set of Y-boundary nodes with zero boundary conditions
  zerobndY1 = bnodeY1( bc_nodeY1 <= 1e-10 );  
  zerobndY2 = bnodeY2( bc_nodeY2 <= 1e-10 );  
  zerobndY3 = bnodeY3( bc_nodeY3 <= 1e-10 );  
  zerobndY  = [zerobndY1; zerobndY2; zerobndY3];
  
% Set of Y-boundary nodes with non-zero boundary conditions  
  nonzerobndY1 = bnodeY1( ~ismember(bnodeY1,zerobndY1) );
  nonzerobndY2 = bnodeY2( ~ismember(bnodeY2,zerobndY2) );
  nonzerobndY3 = bnodeY3( ~ismember(bnodeY3,zerobndY3) );
  nonzerobndY  = [nonzerobndY1; nonzerobndY2; nonzerobndY3];
   
end % end child function


% ------------------------------------------------------------
% Child function
% ------------------------------------------------------------
function [error,bc_nodeY] = interpolated_error_ybc(xy,xyY,bnodesX,bnodeY)
% Impose the boundary conditions on Y-boundary nodes and 
% compute the error

% Get the coordinates of the boundary X-nodes and of the boundary Y-node
  allxbd_X = reshape( xy(bnodesX,1), size(bnodesX,1), 2);
  allybd_X = reshape( xy(bnodesX,2), size(bnodesX,1), 2);
  xybd_Y   = xyY(bnodeY, :);
  
% Compute the boundary conditions for the given nodes
  [bc_firstNodeX]  = specific_bc(allxbd_X(:,1), allybd_X(:,1));
  [bc_secondNodeX] = specific_bc(allxbd_X(:,2), allybd_X(:,2));
  [bc_nodeY]       = specific_bc(xybd_Y(:,1),   xybd_Y(:,2));
 
% Interpolated error
  error = bc_nodeY - 0.5 * (bc_firstNodeX + bc_secondNodeX);

end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [errelem] = get_errelem(evtY,erredges)
%Recovering the element indicators from the 3 edges indicators per element
    
% Get the matrix nel-by-3 edgelem having per each row (element) the 
% correspoding edge-indicator (column)
  evtYvec = reshape(evtY',3*size(evtY,1),1);
  edgelem = erredges(evtYvec);
  edgelem = reshape(edgelem,3,length(edgelem)/3)';
  
% Compute the elementwise estimates by square summing the 3 edge 
% indicators per element
  errelem = sqrt( sum(edgelem.^2, 2) );
          
end % end function 