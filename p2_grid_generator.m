function [p2xy,p2evt,p2bound] = p2_grid_generator(xy,evt,bound)
%P2_GRID_GENERATOR quadratic element grid generator
%
% [p2xy,p2evt,p2bound] = p2_grid_generator(xy,evt,bound)
%
% input:
%               xy      vertex coordinate vector  
%              evt      element vertex mapping matrix
%            bound      boundary vertex vector
% output:
%             p2xy      P2 nodal coordinate vector
%            p2evt      P2 element node mapping matrix
%          p2bound      P2 boundary node vector
%
% Given an input mesh (by xy, evt, and bound), this function returns the p2xy, 
% p2evt, and p2bound data structures containing the mesh xy-coordinates, 
% evt-element map and the vector of boundary nodes associated with a uniform 
% refinement of the mesh. 
%
% Uniform means that a refinement using these data structures would introduce 
% the 3 edge-midpoints per element.
%
% Note that no real refinement is performed here: only coordinates and numbers 
% of the edge midpoints are returned. 
%
% See also UNIFORM_REFINEMENT
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  nvtx = size(xy,1);    % Number of vertices
  nel  = size(evt,1);   % Number of elements 
  
% Total edges (with repetitions)  
  edges = [evt(:,[2,3]); evt(:,[3,1]); evt(:,[1,2])];
  
% Total edges of the triangulation without repetitions
  [edges,~,ic] = unique(sort(edges,2),'rows');    

% Edges per element (i.e., the new grid nodes per element)
  edgeNum = [ic(1:nel), ic(nel+1:2*nel), ic(2*nel+1:3*nel)];
    
% Compute the P2 elements map:
% This is a nel-by-6 matrix containing the P1 evt on the first 3 columns,
% and the associated midpoints on the last 3 columns
  p2nodes = edgeNum + nvtx;
  p2evt   = [evt, p2nodes];
  
% Compute the P2 xy-coordinates (i.e., all edge's midpoints):
% This is a allnodes-by-2 matrix constisting of the P1 xy mesh coordinates
% matrix and the coordinates of new midpoints appended at the end 
  mxy(:,1) = 0.5 * (xy(edges(:,1),1) + xy(edges(:,2),1));
  mxy(:,2) = 0.5 * (xy(edges(:,1),2) + xy(edges(:,2),2));
  p2xy     = [xy; mxy];
  
% Finding the midpoint boundary nodes (i.e., the boundary edges)
  checkBound   = sum( ismember(edges,bound) , 2);
  p2nodesBound = find(checkBound==2) + nvtx;  
  p2bound      = [bound; p2nodesBound];
    
end % end function