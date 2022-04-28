function [nodes] = stochcol_nodes_cc_increasing(levels)
%STOCHCOL_NODES_CC computes the (1D) Clenshaw-Curtis (CC) nodes in [-1,1]
%
% input:
%     levels    vector of levels
%
% output:
%     nodes     cell of CC collocation points
%     
% The function returns the 1D Clenshaw-Curtis interpolation sparse grid
% nodes for a given set of levels 'i' (i>0). 
%
% At level i, the number of collocation points are:
%           
%   mi := 1              if i = 1,
%   mi := 2^(i-1) + 1    if i > 1.
%
% Then, the set of CC nodes are defined as:
%
%   y_1^(i) := 0                         i = 1, 
%   y_j^(i) := -cos( pi*(j-1)/(mi-1) )   i > 1, j = 1,...,mi.
%
%   TIFISS function: LR; 27 September 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if ~isequal(length(levels),nnz(levels))
      % There is some zero in the input sequence of levels
      error('Levels can be only numbers greater or equal than 1');
  end
  
% Allocate memory  
  nodes = cell(length(levels),1);
 
% Assign the only node 0.0 at levels=1    
  [nodes{levels==1}] = deal(0.0);
 
% Define the number of nodes for (levels~=1)
  m = @(i) (2^(i-1) + 1);   

% Define the 1D Clenshaw-Curtis rule (levels~=1)  
  clencurt = @(i) -cos( pi * ( 0:1:(m(i)-1) ) / (m(i)-1) ) ; 
    
% Position of levels not equal to 1  
  levnotone = levels(levels~=1);

% Empty cell yet   
  stillempty = find(levels~=1);
  
  for i = 1:length(levnotone)
      nodes{stillempty(i)} = clencurt( levnotone(i) );
      %
      % Update zeros to 0.0 in the central position
      midpos = ( 1 + m(levnotone(i)) ) / 2;       
      nodes{stillempty(i)}(midpos) = 0.0;
      % This step is needed because -cos(pi/2) ~= 1e-17, which is not
      % recognized as exactly 0
  end

end % end function