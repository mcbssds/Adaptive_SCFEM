function [sigptloc,tigptloc] = subelt_transf(s,t,subelt,subdivPar)
%SUBELT_TRANSF locate 2D Gaussian points on the right reference subelement
%
%  [sigptloc,tigptloc] = subelt_transf(s,t,subelt)
%
%  input:
%             s    x-coordinate 2D Gaussian point in the reference triangle
%             t    y-coordinate 2D Gaussian point in the reference triangle
%        subelt    sub-element
%     subdivPar    uniform or criss-cross subdivision flag
%
%  output:
%        sigptloc  x-coordinate 2D Gaussian point in the right sub-element
%        tigptloc  y-coordinate 2D Gaussian point in the right sub-element
%          
% Last update: 01/02/2017
% -------------------------------------------------------------------------
%    TIFISS function:
% Copyright (c) 2017 Alex Bespalov, Leonardo Rocchi

  if subdivPar == 1 
      % Uniform subdivision
      % -------------------------------------------------------------------
      if     subelt == 1
          % Top sub-element
          x = [  0, 1/2, 0];    % x-coordinates sub-element
          y = [1/2, 1/2, 1];    % y-coordinates sub-element  
      elseif subelt == 2 
          % Left-bottom sub-element
          x = [0, 1/2,   0];   
          y = [0,   0, 1/2]; 
      elseif subelt == 3          
          % right-bottom sub-element  
          x = [1/2, 1, 1/2];
          y = [  0, 0, 1/2];
      else
          % Central sub-element
          x = [1/2, 0, 1/2];
          y = [1/2, 1/2, 0];
      end
      
  else
      % Criss-cross subdivision
      % -------------------------------------------------------------------
      if subelt == 1
          % Top sub-element
          x = [0.0, 0.0, 0.5];
          y = [1.0, 0.5, 0.5];
      elseif subelt == 2
          % Top-left sub-element
          x = [0.5, 0.0, 0.0]; 
          y = [0.5, 0.5, 0.0];
      elseif subelt == 3
          % Bottom-left sub-element
          x = [0.0, 0.5, 0.5]; 
          y = [0.0, 0.0, 0.5];
      else%if subelt == 4
          % Bottom-left sub-element
          x = [0.5, 0.5, 1.0]; 
         y = [0.5, 0.0, 0.0];
      end
      
  end
  
% Affine mapping from the reference element K^ to the right sub-element
  B = [x(2) - x(1),  y(2) - y(1); ...
       x(3) - x(1),  y(3) - y(1)];
  c = [x(1); y(1)];
  F = c + B'*[s;t];

% Assign new points
  sigptloc = F(1);
  tigptloc = F(2);

end  % end function
