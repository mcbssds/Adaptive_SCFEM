function [X] = stochcol_getindexset(w,M,iplot)  
%STOCHCOL_GETINDEXSET computes the set of multi-indices X(w,M)
%
% [X] = stochcol_getindexset(w,M,iplot)
% 
% input:
%        w    level of the sparse grid (integer >= 0)
%        M    number of directions (i.e. no. of (active) parameters noarv)
%    iplot    (optional) switch plot (=1) of the multi-indices in 2D and 3D  
%
% output:    
%        X    set of multi-indices for the sparse grid 
%
% The function returns the multindex set X(w,M) for the associated sparse grid 
% defined by 
%
%   X(w,M) := { ii=(i_1,...,i_M) : sum_{k=1}^{M}(ii_k - 1) <= w }
%
% which corresponds to eq. (2.7a) in the main reference:
% * [NTW08] Nobile, Tempone, Webster, 'A sparse grid stochastic collocation 
% * method for partial differential equations with random input data', SIAM
% * J. Numer. Anal, 46(5), pp. 2309-2345, 2008. 
%
% The multindex set is computed using the function MULTIIDX_GEN for a
% general "rule" from the:
% * Sparse grids Matlab kit *
% * by Back, Nobile, Tamellini, and Tempone *
% * available at: https://csqi.epfl.ch/page-107231-en.html *
% See also the function FAST_TD_SET from the kit which directly computes X(w,M) 
% for CC-rule.
%
% Calling as:
% - [X] = stochcol_getindexset(w,M)    (no plot)
% - [X] = stochcol_getindexset(w,M,1)  (plots the sparse grid (if M=2/3))
%
% Function(s) called: multiidx_gen
%
%   TIFISS function: LR; 27 September 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if w < 0 && M < 1
      error('Wrong input arguments!');
  elseif w < 0
      error('First argument has to be greater of equal than 0!');
  elseif M < 1
      error('Second argument has to be greater or equal than 1!'); 
  end
  
  if nargin < 3 
      iplot = 0;
  end

% Calling STOCHCOL_MULTIIDX_GEN, a copy of original function from "Sparse grids Matlab kit"
% The default rule is rule=@(i) sum(i-1), see definition of X(w,X) above.
% This is also the rule used for the CC-rule (see also FAST_TD_SET);
  rule=@(I) sum(I-1);  
  X = stochcol_multiidx_gen(M,rule,w,1);
  
  if iplot%==1  
      if M==2
          figure;
          scatter(X(:,1),X(:,2),100,'o','filled');
          title(['Multindex set X(',num2str(w),',',num2str(M),')  =  ( i  :  i_1 + i_2  <=  ',num2str(w+M),')']);
          grid on; axis square;
          xlabel('i_1'); ylabel('i_2');
          axis([ 0 max(X(:,1))+1 0 max(X(:,1))+1 ]);
          set(gca,'XTick',0:max(X(:,1))+1,'YTick',0:max(X(:,1))+1,'Fontsize',14);
      
      elseif M==3
          figure;
          scatter3(X(:,1),X(:,2),X(:,3),100,'o','filled');
          title(['Multindex set X(',num2str(w),',',num2str(M),')  =  ( i  :  i_1 + i_2 + i_3 <=  ',num2str(w+M),')']);
          grid on; axis square;
          view(50,20);
          xlabel('i_1'); ylabel('i_2'); zlabel('i_3');
          axis([ 1 max(X(:,1))+1 1 max(X(:,1))+1 1 max(X(:,1))+1]);
          set(gca,'XTick',0:max(X(:,1))+1,'YTick',0:max(X(:,1))+1,'ZTick',0:max(X(:,1))+1,'Fontsize',14);
      end
  end
    
end % end function