function [Mset] = marking_strategy_fa(indicators,strategy,threshold)
%MARKING_STRATEGY_FA marking strategy for a vector of error indicators
%
% [Mset] = marking_strategy_fa(indicators,strategy,threshold)
%
% input:
%       indicators   vector of indicators (can be associated with elements/edges...)
%         strategy   marking strategy
%        threshold   threshold marking parameter
% output:
%             Mset   set of marked elements associated with indicators
%
% Strategies implemented are:
% 1 - maximum;
% 2 - Doerfler; 
% 3 - maximum with percentage;
% 4 - Doerfler with percentage.
% Option 3 (resp. 4) first marks a certain percentage (epsilon) of indicators 
% and then applies the associated maximum (resp. Doerfler) strategy to the 
% remanining (i.e., not marked) indicators (more indicators are marked overall).
% Note that these strategies are sometimes more efficient depending on the problem.
% 
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if nargin < 3
      % Default marking parameter
      threshold = 0.5;
      if nargin < 2
          % Default marking strategy: equilibration strategy
          strategy = 2;
          if nargin < 1
              error('Insufficient input arguments!');
          end
      end
  end
        
% Number of elements
  lenind = length(indicators); 

% Choose strategy
  if strategy == 1
      % --------------------------------------------------------------
      % Maximum strategy
      % --------------------------------------------------------------
      max_err = max(indicators);
      Mset = find(indicators >= threshold * max_err);
      %Mset = sort(Mset);
   
  elseif strategy == 2
      % --------------------------------------------------------------
      % Doerfler (equilibration) strategy
      % --------------------------------------------------------------
      if isequal(threshold,0)
          % Empty set
          Mset = [];
      elseif isequal(threshold,1)
          % Full set
          Mset = (1:lenind)';
      else
          % Sorting indicators (from larger to smaller)
          [sortinds,idx] = sort(indicators,'descend');
          % Cumulative sum (cumsum) of indicators^2
          cumerrs = cumsum( sortinds.^2 );            
          % Set of marked elements/edges/indices
          Mset = idx( 1 : find(cumerrs >= threshold * sum(indicators.^2),1) );
      end
                      
  elseif strategy == 3 
      % --------------------------------------------------------------
      % Maximum strategy with modification percentage
      % --------------------------------------------------------------
      epsilon = 0.1;
      thrs    = floor(lenind*epsilon);
      [~,I]   = sort(indicators,'descend');
      % First set M1 
      M1 = I(1:thrs);
      %
      % Now apply the maximum strategy only for the indicators which are not 
      % in M1: this is equivalent to consider indicators(M1)=0.
      indicators(M1) = 0.0;
      max_err = max(indicators);
      M2 = find(indicators >= threshold * max_err);
      %
      % Union with the first set M1 
      Mset = [M1;M2];  %sort([M1;M2]);
   
  elseif strategy == 4
      % --------------------------------------------------------------
      % Doerfler (equilibration) strategy with modification percentage
      % --------------------------------------------------------------
      if isequal(threshold,0)
          % Empty set
          Mset = [];
      elseif isequal(threshold,1)
          % Full set
          Mset = (1:lenind)';
      else
          epsilon = 0.2;
          thrs = floor(lenind*epsilon);
          % Sorting the estimators
          [values,idx] = sort(indicators,'descend');
          % First set of marked indicators
          M1 = idx(1:thrs);
          %
          % Now apply the Doerfler strategy to the remaining indicators
          %          
          % Second set where to mark the others
          indic2 = values(thrs+1:end);
          %
          % Sorting indicators (from larger to smaller)
          [sortinds2,idx2] = sort(indic2,'descend');
          % Cumulative sum (cumsum) of indicators^2
          cumerrs = cumsum( sortinds2.^2 );            
          % Set of marked elements/edges/indices
          M2 = idx2( 1 : find(cumerrs >= threshold * sum(indic2.^2),1) );
          %
          % Union with the first set
          Mset = [M1;M2];  %sort([M1;M2]);    
      end
  
  end % end if

end % end function