      function [xi,dxids,dxidt] = tshape(s,t)
%TSHAPE  evaluates linear shape functions
%   [xi,dxids,dxidt] = tshape(s,t);
%   input
%          s         first triangle coordinate   
%          t         second triangle coordinate  
%  
%   output
%          xi        shape function
%          dxids     s derivative of xi
%          dxidt     t derivative of xi
%
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
      one = 1.0e0; zero=0.0e0;
      xi(1) = one-s-t;
      xi(2) = s;
	  xi(3) = t;
      dxids(1) = -one;
      dxids(2) = one;
      dxids(3) = zero;
      dxidt(1) = -one;
      dxidt(2) = zero;
      dxidt(3) = one;
      return