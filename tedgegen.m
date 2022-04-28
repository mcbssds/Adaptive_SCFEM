function [eex,tve,els] = tedgegen(xy,evt)
%TEDGEGEN edge information for flux jump computation
%
%   [eex,tve,els] = tedgegen(xy,evt)
%
%   input:
%             xy     nodal coordinate vector 
%            evt     element mapping matrix
%
%   output:
%            eex     element-connectivity matrix
%            tve     edge-location matrix
%            els     elementwise edge lengths
%
%   TIFISS function: LR; 05 October 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
 
  nel  = size(evt,1);   % Number of elements
  nvtx = size(xy,1);    % Number of vertices
  x = xy(:,1);
  y = xy(:,2);
   
% Recover local coordinates
  xl_v = zeros(nel,3);
  yl_v = zeros(nel,3);
  for ivtx = 1:3
      xl_v(:,ivtx) = x(evt(:,ivtx));
      yl_v(:,ivtx) = y(evt(:,ivtx));
  end

% Compute edge lengths
% ------------------------------------------------------------------
% First edge
  hx_v = xl_v(:,3) - xl_v(:,2); 
  hy_v = yl_v(:,3) - yl_v(:,2);
  els(:,1) = sqrt(hx_v.^2 + hy_v.^2);
 
% Second edge
  hx_v = xl_v(:,1) - xl_v(:,3); 
  hy_v = yl_v(:,1) - yl_v(:,3);
  els(:,2) = sqrt(hx_v.^2 + hy_v.^2);

% Third edge
  hx_v = xl_v(:,2) - xl_v(:,1); 
  hy_v = yl_v(:,2) - yl_v(:,1);
  els(:,3) = sqrt(hx_v.^2 + hy_v.^2);

  %fprintf('mesh validation step\n')
% Validate the mesh
  detel = (xl_v(:,2) .* yl_v(:,3) - yl_v(:,2) .* xl_v(:,3)) ...
        - (xl_v(:,1) .* yl_v(:,3) - yl_v(:,1) .* xl_v(:,3)) ...
        + (xl_v(:,1) .* yl_v(:,2) - yl_v(:,1) .* xl_v(:,2));
  if any(detel <= 0)
      error('Oops...mesh issue: look at element %d',find(detel<=0));
  end

% Create the adjacency matrix of pairs of elements sharing the same edge  
% Ex:      (i,j) = 3   and   (j,i) = 5
% It means that elements 3 and 5 share the edge i-j (or j-i). the matrix
% has "symmetric" positions since edge i->j is, for example, the 1st edge
% of element 3, and j->i is the 3rd edge of element 5.
  adj = sparse(nvtx,nvtx);
  adj = adj + sparse(evt(:,2), evt(:,3), 1:nel, nvtx,nvtx);   % All the first edges
  adj = adj + sparse(evt(:,3), evt(:,1), 1:nel, nvtx,nvtx);   % All the second edges
  adj = adj + sparse(evt(:,1), evt(:,2), 1:nel, nvtx,nvtx);   % All the third edges
 
% Create the element-connectivity matrix eex (nel-by-3)
% Ex 1:    row i-th    ->    14    6    9
% It means that the neighbours of element i-th are elements 14, 6, and 9, 
% and they lie respectively on the 1st, 2nd, and 3rd edge 
% (they are "sorted" with respect to edges)
% Ex 2:    row i-th    ->    14    0    9
% It means that element i-th is a boundary element and the boundary edge
% is the 2nd one.
  idx1  = {evt(:,3), evt(:,2)};
  ind2  = {evt(:,1), evt(:,3)};
  ind3  = {evt(:,2), evt(:,1)};
  elem1 = full( adj( sub2ind(size(adj),idx1{:}) ) );
  elem2 = full( adj( sub2ind(size(adj),ind2{:}) ) );    
  elem3 = full( adj( sub2ind(size(adj),ind3{:}) ) );
  eex   = [elem1, elem2, elem3];
      
% Replace zeros with the same boundary elements of the corresponding rows:
% this is needed otherwise we cannot create the next matrix edloc
  [ii,jj] = find(eex == 0);
  idx = {ii, jj};
  eex( sub2ind(size(eex),idx{:}) ) = ii;
  
% Create the adjacency edge-location matrix edloc (nel-by-nel sparse)
% -------------------------------------------------------------------------
% Ex:   row i-th  ->   0   0   (i,p)=3   0   0   (i,q)=1   0   (i,r)=2
% It means that the neighbours of the i-th element are elements p,q,r and 
% - element p lies on the 3rd edge of i;
% - element q lies on the 1st edge of i;
% - element r lies on the 2nd edge of i.
% -------------------------------------------------------------------------
% If we look at edloc column-wise...
% Ex:   col j-th  ->   0   0   (p,j)=1   0   0   (q,j)=3   0   (r,j)=2
% It means that the neighbours of the j-th element are elements p,q,r and
% - element j lies on the 1st edge of element p;
% - element j lies on the 3rd edge of element q;
% - element j lies on the 2nd edge of element r.
% Looking at edloc column-wise we can construct the tve matrix:
  edloc = sparse(nel,nel);
  for j = 1:3
      edloc = edloc + sparse((1:nel)', eex(:,j), j*ones(nel,1), nel,nel);
  end
    
% The boundary elements have the boundary edges along the main diagonal: 
% replace them with zeros now
  idx = {(1:nel)', (1:nel)'};
  edloc( sub2ind(size(edloc), idx{:}) ) = 0;

% Create the edge-location matrix tve (nel-by-3)
% Ex 1:   row i-th    ->    1    3    2
% It means that i-th element lies respectively on the 1st, 3rd, and 2nd 
% edge of its corresponding neighbours, which are stored in eex(i,:)
% Ex 2:   row i-th    ->    0    2    1
% It means that i-th element is a boundary element (in particular, its 
% boundary edge is the 1st one). Then, it has to be read as in Ex 1  
  ii = reshape(eex', 3*nel, 1);               % row-indices
  jj = reshape(repmat(1:nel,3,1), 3*nel, 1);  % col-indices
  idx = {ii, jj};
  tve = full( edloc( sub2ind(size(edloc), idx{:}) ) );
  tve = reshape(tve,3,nel)';

end % end function
