function [elerr, ederr, errest] = diffpost_p1_with_p1_linsys(evt, xy, ...
    eboundt, xgal, evtY, xyY, boundY, Ybasis, subdivPar, coeff_fun, ...
    rhs_fun)
%DIFFPOST_P1_WITH_P1_LINSYS computes hierarchical eY estimator solving the (fully) assembled error problem
%
%[elerr, ederr, errest] = diffpost_p1_with_p1_linsys(evt, xy, ...
%    eboundt, xgal, evtY, xyY, boundY, Ybasis, subdivPar, coeff_fun, ...
%    rhs_fun)
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
%            elerr    vector of hierarchical element indicators
%            ederr    vector of hierarchical edge indicators
%           errest    global hierarchical error estimate
%
% The function solves the discrete problem for the eY hierarchical estimator:
%
%   B(eY,v) = F(v) - B(xgal,v),     for all v \in Y,        (1)
%
% NOTE that given the midpoint-solutions (edge solutions) to problem (1) per 
% each node, the element indicators are computed by definition of local 
% B-norms restricted to single elements.
%                    
% Function(s) called: triangular_gausspoints
%                     tderiv
%                     tgauss
%                     subelt_transf
%                     specific_bc
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  nvtx = size(xy,1);      % Number of vertices  (i.e., number of X-basis functions) 
  nyb  = size(Ybasis,1);  % Number of midpoints (i.e., number of Y-basis functions (boudary ones included))

% Extract elements and associated edges for Y-basis functions
  elems = Ybasis(:,[1,2]);   
  edges = Ybasis(:,[3,4]);  
 
% -----------------------------------------------------------------------------                 
% STEP 1: elementwise contributions from uniform red/bisec3 refinement type
% -----------------------------------------------------------------------------      
  [ae,bde,~,fde] = diffpost_p1_detcontrib(xy,evt,subdivPar,coeff_fun, rhs_fun);
  
% -----------------------------------------------------------------------------
% STEP 2: assembling the LHS matrix
% -----------------------------------------------------------------------------
  A = sparse(nyb,nyb);
  for i = 1:3
      % Indices for Y-basis 
      idxY1 = evtY(:,i);
      for j = 1:3
          % Indices for Y-basis
          idxY2 = evtY(:,j);
          % Assembling
          A = A + sparse(idxY1, idxY2, ae(:,i,j), nyb, nyb);
      end
  end  
       
% -----------------------------------------------------------------------------                       
% STEP 3: assembling F(v) for all basis v in Y
% -----------------------------------------------------------------------------
% Indices: elements, Y-basis
  ind1 = {elems(:,1), edges(:,1)};
  ind2 = {elems(:,2), edges(:,2)};
% Summing contributions
  ff = fde( sub2ind(size(fde) , ind1{:}) ) + fde( sub2ind(size(fde) , ind2{:}) );
% For a Y-basis on the boundary, we have to halve the contributions, since in the 
% above line we doubled the contributions coming from the boundary elements
  ff(boundY) = ff(boundY) ./ 2; 
  
% -----------------------------------------------------------------------------                                
% STEP 4: Assembling B(ugal,v) for all basis v in Y
% -----------------------------------------------------------------------------  
% The term B(ugal,v) will be a vector nyb-by-1     
  Km = sparse(nyb,nvtx);
  for kYb = 1:3
      indY = evtY(:,kYb);	 
      for kXb = 1:3
          indX = evt(:,kXb);
          Km = Km  + sparse(indY, indX, bde(:, kYb, kXb), nyb, nvtx);
      end
  end  
  Bx = Km * xgal;
      
% Global RHS of the error problem: F(v) - B(ugal,v)  
  rhs = ff - Bx;
   
% -----------------------------------------------------------------------------
% STEP 5: imposing interpolated error as Dirichlet bcs on Y-boundary nodes
% -----------------------------------------------------------------------------
  [rhs] = diffpost_nonzerobc(evt,evtY,xy,xyY,boundY,eboundt,A,rhs);

% -----------------------------------------------------------------------------
% STEP 6: solving the associated assembled system for interior nodes
% -----------------------------------------------------------------------------
% The final ederr will contain the error solutions at each midpoint.
% NOTE that they are nothing but that the edge (or midpoint) indicators.
  ederr = zeros(nyb,1);

% Set of interior midpoints 
  interiorY = 1:boundY(1)-1;
   
% Solving the linear systems for interior positions
  ederr(interiorY) = A(interiorY,interiorY) \ rhs(interiorY);

% Update the error values in boundary-midpoints position
  ederr(boundY) = rhs(boundY);

% -----------------------------------------------------------------------------
% STEP 7: recovering the element indicators from edge indicators
% -----------------------------------------------------------------------------  
  [elerr] = get_errelem(ae,evtY,ederr);
  
% -----------------------------------------------------------------------------  
% STEP 8: global estimates ||eY||_B0^2 = B0(eY,eY)
% -----------------------------------------------------------------------------  
% Compute ||eY||_B0^2 = B0(eY,eY)
  errest_sq = ederr' * rhs;
  
% Compute ||eY||_B0 = B0(eY,eY)^(1/2)
  errest = full( sqrt( errest_sq ) );
  
end % end function



% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [ade, bde, bdem, fde] = diffpost_p1_detcontrib(xy, evt, subdivPar, coeff_fun, rhs_fun)
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
%   B(eY,v) = F(v) - B(xgal,v) for all v \in Y
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
% STEP 3: integrals in B(eY,v) for all basis v in Y
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
          %
          % Evaluate derivatives
          [~,invjac_v,~,dphidx_v,dphidy_v] = tderiv(sigpt,tigpt,xl_m,yl_m);
          %
          % Evaluate variable diffusion coefficients
          [coeff] = tgauss(sigpt, tigpt, xl_m, yl_m, coeff_fun);
          %
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
function [rhs] = diffpost_nonzerobc(evt,evtY,xy,xyY,boundY,eboundt,K,rhs)
%Imposes Dirichlet boundary conditions at boundary midpoints
  
  interiorY = 1:boundY(1)-1;

% -----------------------------------------------------------------------------
% Extract boundary elements with boundary edges respectively = 1, 2, and 3
% ----------------------------------------------------------------------------- 
  beled1 = eboundt( eboundt(:,2)==1 , 1);     %nbel1 = length(beled1);
  beled2 = eboundt( eboundt(:,2)==2 , 1);     %nbel2 = length(beled2);
  beled3 = eboundt( eboundt(:,2)==3 , 1);     %nbel3 = length(beled3);
   
% -----------------------------------------------------------------------------
% Extract boundary nodes/midpoints and compute the error
% -----------------------------------------------------------------------------
  node1    = evt ( beled1, [2 3]);   % nodes of the 1st boundary edges (columns 2,3 of evt)
  midp1    = evtY( beled1, 1);       % midpoints on 1st boundary edges (column 1 of evtY)
  [error1] = interperror_bc(xy,xyY,node1,midp1);
  
  node2    = evt ( beled2, [3 1]);   % nodes of the 2nd boundary edges (columns 3,1 of evt)
  midp2    = evtY( beled2, 2);       % midpoints on 2nd boundary edges (column 2 of evtY)
  [error2] = interperror_bc(xy,xyY,node2,midp2);
  
  node3    = evt ( beled3, [1 2]);   % nodes of the 3rd boundary edges (columns 1,2 of evt)
  midp3    = evtY( beled3, 3);       % midpoints on 3rd boundary edges (column 3 of evtY)
  [error3] = interperror_bc(xy,xyY,node3,midp3);
  
% -----------------------------------------------------------------------------
% Vector of interpolated error on boundary midpoints
% -----------------------------------------------------------------------------
  miderrors = sortrows([midp1, error1; midp2, error2; midp3, error3],1);
  bcY = miderrors(:,2); 
  
% -----------------------------------------------------------------------------
% Update the Rhs in the interior-midpoint positions
% -----------------------------------------------------------------------------
% We have to do the following:
%
% rhs(interiorY) = rhs(interiorY) - K(interiorY,boundY) * bcY
%
% Note that this updated rhs for interior positions is a nIntY-by-1 vector;
  rhs(interiorY)  = rhs(interiorY) - K(interiorY,boundY) * bcY;

% Update values in boundary-midpoints position with the imposed error
  rhs(boundY)  = bcY;
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [error] = interperror_bc(xy,xyY,bnodesX,bnodeY)
%Computes the error committed at boundary midpoints

% Coordinates of the boundary nodes and of boundary edge's midpoints
  allxbd_X = reshape( xy(bnodesX,1), size(bnodesX,1), 2);
  allybd_X = reshape( xy(bnodesX,2), size(bnodesX,1), 2);
  xybd_Y   = xyY(bnodeY, :);
  
% Compute boundary conditions for the given nodes
  [bc_firstNodeX]  = specific_bc( allxbd_X(:,1), allybd_X(:,1) );
  [bc_secondNodeX] = specific_bc( allxbd_X(:,2), allybd_X(:,2) );
  [bc_nodeY]       = specific_bc( xybd_Y(:,1),   xybd_Y(:,2)   );

% Interpolated error
  error = bc_nodeY - 0.5 * (bc_firstNodeX + bc_secondNodeX);
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [elerr] = get_errelem(ade,evtY,edgevec)
%Recovers elementwise B_0^2 energy norms

  nel      = size(evtY,1);
  errloc   = zeros(nel,3);
  elerr_sq = zeros(nel,1);
     
% Rearranging spatial nodal values elementwise
  for imid = 1:3
      errloc(:,imid) = edgevec(evtY(:,imid));
  end 
     
% Computing elementwise B_0^2 energy norms
  for imid = 1:3          
      for jmid = 1:3
          elerr_sq(:) = elerr_sq(:) + errloc(:,imid) .* ade(:,imid,jmid) .* errloc(:,jmid);
      end
  end
  
% B_0 energy norms
  elerr = sqrt( elerr_sq );
  
end % end child function 