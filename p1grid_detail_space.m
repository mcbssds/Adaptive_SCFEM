function [evtY,xyY,boundY,Ybasis] = p1grid_detail_space(xy,evt)
%P1GRID_DETAIL_SPACE linear detail space Y grid-generator
%
% [evtY,xyY,boundY,Ybasis] = p1grid_detail_space(xy,evt)
%
% input:
%            xy     vertex coordinate vector  
%           evt     element mapping matrix
%
% output:
%          evtY     element mapping matrix for midpoints
%           xyY     vertex coordinate vector for midpoints
%        boundY     boundary midpoints vector
%        Ybasis     Y-basis element-positions matrix
%
% The function creates data structures corresponding to the detail space Y 
% spanned by piecewise linear basis associated with edge midpoints of the 
% current triangulation. In particular:
% - like evt, evtY is a nel-by-3 matrix containing midpoints per element;
% - Ybasis is #allmidpoints-by-4 matrix containing the elements sharing each 
%   midpoint (neighbours) as well as the corresponding edge-position;
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi
                          
  nel  = size(evt,1);   % Number of elements
  nvtx = size(xy,1);    % Number of vertices

% -------------------------------------------------------------------------
% STEP 1: create the Ybasis matrix (see description below)
% -------------------------------------------------------------------------
  
% Create the adjacency matrix of pairs of elements sharing the same edge  
% Ex:    (i,j) = 3   and   (j,i) = 5
% It means that elements 3 and 5 share the edge i-j (or j-i). The matrix
% has "symmetric" positions since edge i->j is, for example, the 1st edge
% of element 3, and j->i is the 3rd edge of element 5.
  adj = sparse(nvtx,nvtx);
  adj = adj + sparse(evt(:,2), evt(:,3), 1:nel, nvtx,nvtx);   % All the first edges
  adj = adj + sparse(evt(:,3), evt(:,1), 1:nel, nvtx,nvtx);   % All the second edges
  adj = adj + sparse(evt(:,1), evt(:,2), 1:nel, nvtx,nvtx);   % All the third edges
 
% Create the element-connectivity-matrix eex (nel-by-3)
% Ex:    i-th row    ->    14    6    9
% It means that the neighbours of the i-th element are elements 14, 6, and 9
% and they lie on the 1st, 2nd, and 3rd edge, respectively (they are "sorted" 
% with respect to edges;
% Ex:    i-th row    ->    14    0    9
% It means that the i-th element is a boundary element and the boundary edge
% is the 2nd one.
  idx_one    = {evt(:,3), evt(:,2)};
  ind_two    = {evt(:,1), evt(:,3)};
  ind_three  = {evt(:,2), evt(:,1)};
  elem_one   = full( adj( sub2ind(size(adj),idx_one{:}) ) );  
  elem_two   = full( adj( sub2ind(size(adj),ind_two{:}) ) );   
  elem_three = full( adj( sub2ind(size(adj),ind_three{:}) ) );
  eex        = [elem_one, elem_two, elem_three];
  
% Replace zeros with the same boundary elements of the corresponding rows:
% this is needed in order to create the next matrix 'edloc'
  [ii,jj] = find(eex == 0);
  idx     = {ii, jj};
  eex( sub2ind(size(eex),idx{:}) ) = ii;
  
% Create the adjacency edge-location matrix edloc (nel-by-nel sparse)
% Ex:    i-th row  ->   0   0   (i,p)=3   0   0   (i,q)=1   0   (i,r)=2
% It means that the neighbours of the i-th element are elements p,q,r and 
% - element p lies on the 3rd edge of i;
% - element q lies on the 1st edge of i;
% - element r lies on the 2nd edge of i;
% If we look at edloc column-wise:
% Ex:    j-th col ->   0   0   (p,j)=1   0   0   (q,j)=3   0   (r,j)=2
% It means that the neighbours of the j-th element are elements p,q,r and
% - element j lies on the 1st edge of element p;
% - element j lies on the 3rd edge of element q;
% - element j lies on the 2nd edge of element r.
  edloc = sparse(nel,nel);
  for j = 1:3
      edloc = edloc + sparse((1:nel)', eex(:,j), j*ones(nel,1), nel,nel);
  end
% Boundary elements have the boundary edges along the main diagonal
  
% Take the upper part of edloc excluding the main diagonal (boundary edges)
  uppPart = triu(edloc,1);
  [elem1,elem2,pos1] = find(uppPart);
  
% Find the position with respect to elements in elem2 by taking the entries  
% in the lower symmetric positions of edloc
  idx  = {elem2, elem1};
  pos2 = full( edloc( sub2ind(size(edloc),idx{:}) ) ); 
  
% Final Ybasis matrix (#internaledges-by-4)
% Ex:    i-th row    ->    (i)    7    6    3    2
% It means that the i-th midpoint (Y-basis) will be shared by elements 7 and 6,
% and it will lie on the 3rd edge of element 7 and on the 2nd edge of element 6. 
% NOTE that at this stage Ybasis does not include the boundary Y-basis 
% functions (midpoints); it will updated later.
  Ybasis = [elem1, elem2, pos1, pos2];
   
% -------------------------------------------------------------------------
% STEP 2: create the uniform refinement mapping matrix evtY
% -------------------------------------------------------------------------
% Ex:    i-th row    ->     8    3    14
% It means that the Y-basis functions (midpoints) of the i-th element are 
% 8, 3, 14 and the lie on the 1st, 2nd, 3rd edge of the element, respectively;
% it is similar to evt but containing midpoints instead of vertices.
% Ex:    i-th row    ->     8    0    14
% It means that the 2nd Y-basis function (midpoint) will be on the boundary
  nYb  = size(Ybasis,1); 
  evtY = zeros(nel,3);
  idx1 = {elem1, pos1};
  idx2 = {elem2, pos2};
  evtY( sub2ind(size(evtY), idx1{:}) ) = (1:nYb);
  evtY( sub2ind(size(evtY), idx2{:}) ) = (1:nYb); 
  
% -------------------------------------------------------------------------
% STEP 3: update the Ybasis and evtY matrices
% -------------------------------------------------------------------------
% We want to update Ybasis and evtY including the boundary Y-basis
% functions (midpoints):
% - in evtY, they have to be inserted where currently there are zeros: new 
%   basis are inserted starting from the last internal Y-basis function;
% - then, we update Ybasis by simply appending the new rows for boundary 
%   basis to the current YInfo matrix: it will be then #alledges-by-4.

% Number of internal Y-basis 
  nIntBasY = size(Ybasis,1);   
  
% Boundary elements
  [bel,beledge] = find(evtY == 0);
  idx = {bel, beledge};

% Update evtY
  boundY = nIntBasY+1 : nIntBasY+length(bel);
  evtY( sub2ind(size(evtY),idx{:}) ) = boundY; 
 
% Create new rows for the Ybasis-matrix corresponding to the boundary basis
  newrows = [bel, bel, beledge, beledge]; 

% Append rows  
  Ybasis = [Ybasis; newrows];
  
% -------------------------------------------------------------------------
% STEP 4: create the xyY matrix for the coordinates of the midpoints
% -------------------------------------------------------------------------
  x = xy(:,1);
  y = xy(:,2);
  
% Extract Y-basis according to their edge-position (1,2 or 3) w.r.t the 
% elements in the first column of YInfo
  [elEdge1,~] = find(Ybasis(:,3) == 1);
  [elEdge2,~] = find(Ybasis(:,3) == 2);
  [elEdge3,~] = find(Ybasis(:,3) == 3);

% Recovering the 1st, 2nd, 3rd edges w.r.t. to the above Y-basis
  firstEd  = evt( Ybasis(elEdge1,1), [2,3] );   % first edges
  secondEd = evt( Ybasis(elEdge2,1), [3,1] );   % second edges
  thirdEd  = evt( Ybasis(elEdge3,1), [1,2] );   % third edges

% Check that size(frstEd,1) + size(scndEdEd,1) + size(thirdEd,1) == size(YInfo,1)
  if ~isequal( size(firstEd,1) + size(secondEd,1) + size(thirdEd,1), size(Ybasis,1))
      error('Something wrong with the xyY matrix...');
  end
 
% Compute midpoint's coordinates  
  coordEdge1(:,1) = 0.5 * ( x(firstEd(:,1))  + x(firstEd(:,2)) );    
  coordEdge1(:,2) = 0.5 * ( y(firstEd(:,1))  + y(firstEd(:,2)) );
  
  coordEdge2(:,1) = 0.5 * ( x(secondEd(:,1)) + x(secondEd(:,2)) );    
  coordEdge2(:,2) = 0.5 * ( y(secondEd(:,1)) + y(secondEd(:,2)) );
  
  coordEdge3(:,1) = 0.5 * ( x(thirdEd(:,1))  + x(thirdEd(:,2)) );    
  coordEdge3(:,2) = 0.5 * ( y(thirdEd(:,1))  + y(thirdEd(:,2)) );
  
% Create the xy coordinate matrix for the triangulation associated 
% with the detail space Y  
  xyY(elEdge1,:) = coordEdge1;
  xyY(elEdge2,:) = coordEdge2;
  xyY(elEdge3,:) = coordEdge3;

end % end function