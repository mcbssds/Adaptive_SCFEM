function [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mset,imark)
%GET_ALL_MARKED_ELEM recovers set of marked elements/edges (avoiding hanging nodes)
%
% [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mset,imark)
%
% input:
%      Ybasis   Y-basis element-positions matrix
%        evtY   element mapping matrix for midpoints
%        Mset   vector of marked objects, either elements or edges
%       imark   1 or 2 if Mset is a vector of elements or edges, respectively
%
% output:
%       MMele   vector of overall marked elements (hanging nodes avoided)
%      MMedge   vector of overall marked edges    (hanging nodes avoided)
%
% Given in input either the set of marked elements or edges, the function 
% returns both the sets of *overall* marked elements and edges such that 
% the refinement does not produce any hanging node. 
%
% NOTE that *no mesh* refinement is performed here!
%
% NOTE that the function has to be called as:
% - [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Mele,1), 
%   if Mele is the set of marked elements;
%
% - [MMele,MMedge] = get_all_marked_elem(Ybasis,evtY,Medge,2),
%   if Medge is the set of marked edges;
%
% The function works according to the following procedure:
%
% (1) Recover the elements (neighbours) sharing (the midpoints of) the longest 
%     edges of the marked elements Mele: -> Ybasis( evtY(Mele,2) ,[1,2]);
%
% (2) If this set is exactly Mele, it means that no hanging nodes will be
%     produced and no further elements need to be marked. In particular, 
%     this happens if Mele consists of elements sharing only the longest edge;
%
% (3) Otherwise, there is some edge's midpoint belonging to an element
%     which has not been marked. Hence, that midpoint will be an hanging node. 
%     In this case, we add the associated element to the vector of marked 
%     elements.
%
% The procedure continues until condition (2) is satisfied.
%
%   TIFISS function: LR; 04 July 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if ~ismember(imark,[1,2])
      error('Fourth argument either equal to 1 or 2!'); 
  end
  
  if imark==1
      % The input argument is a set of marked elements:
      % - assign Mset to Mele;
      % - set Medge=empty;
      Mele  = Mset;
      Medge = [];      
  else%imark==2
      % The input argument is a set of marked edges:
      % - assign Mset to Medge; 
      % - recover the associated set of "marked elements" which share the marked edges
      Medge = Mset;
      Mele  = Ybasis(Medge,[1 2]);
      Mele  = unique( reshape(Mele,2*size(Mele,1),1) );
  end
  
% Marked elements are saved in a cell: this improve efficiency when new marked 
% elements are added during the loop below since the vector is directly replaced 
% by the new one; in particular, neither preallocation nor 'sparse' is required 
  MMele    = cell(1,1);
  MMele{1} = Mele; 
  
% Updating set  
  while true 
      %
      % Midpoints of longest edges of the current set of marked elements  
      midPointsLe = evtY(MMele{1},2);  
      %
      % Neighbouring-elements sharing midPointsLe
      neighElem = Ybasis( midPointsLe,[1,2] );
      %
      % New marked elements: they are those elements that do not belong to Mele
      % but some of their edges has to be bisected in order to keep conformity
      newMarkdElem = unique( neighElem( ~ismember(neighElem , MMele{1}) ) );
      %          
      % Check if the conformity is recovered
      if isempty(newMarkdElem)
          break; 
      end
      %   
      % Update the current set of marked elements
      MMele{1} = [MMele{1}; newMarkdElem];
  end
  
% Overall set of marked elements 
  MMele = MMele{1};

% Overall set of marked edges (avoiding repetitions)
  MMedge = unique( [Medge; evtY(MMele,2)],'stable' );
  
end % end function
