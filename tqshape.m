function [psi,dpsids,dpsidt] = tqshape(s,t)
%TQSHAPE  evaluates super-quadratic shape functions
% [psi,dpsids,dpsidt] = tqshape(s,t);
%   input
%          s         first triangle coordinate   
%          t         second triangle coordinate  
%  
%   output
%          psi        shape function
%          dpsids     s derivative of psi
%          dpsidt     t derivative of psi
%
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Qifeng Liao
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
%
 	  psi(1) = (2*xi(1)-1)*xi(1);
	  psi(2) = (2*xi(2)-1)*xi(2);
	  psi(3) = (2*xi(3)-1)*xi(3);
	  psi(4) = 4*xi(2)*xi(3);
	  psi(5) = 4*xi(1)*xi(3);
	  psi(6) = 4*xi(1)*xi(2);
      psi(7) = 27*xi(1)*xi(2)*xi(3);
	  dpsids(1) = (4*xi(1)-1)*dxids(1);
	  dpsids(2) = (4*xi(2)-1)*dxids(2);
	  dpsids(3) = (4*xi(3)-1)*dxids(3);
	  dpsids(4) = 4*(xi(2)*dxids(3) + xi(3)*dxids(2)); 
	  dpsids(5) = 4*(xi(3)*dxids(1) + xi(1)*dxids(3));
	  dpsids(6) = 4*(xi(1)*dxids(2) + xi(2)*dxids(1));
	  dpsids(7) = 27*(xi(2)*xi(3)*dxids(1) + xi(3)*xi(1)*dxids(2) +...
                      xi(1)*xi(2)*dxids(3)); 
	  dpsidt(1) = (4*xi(1)-1)*dxidt(1);
	  dpsidt(2) = (4*xi(2)-1)*dxidt(2);
	  dpsidt(3) = (4*xi(3)-1)*dxidt(3);
	  dpsidt(4) = 4*(xi(2)*dxidt(3) + xi(3)*dxidt(2)); 
	  dpsidt(5) = 4*(xi(3)*dxidt(1) + xi(1)*dxidt(3));
	  dpsidt(6) = 4*(xi(1)*dxidt(2) + xi(2)*dxidt(1));
	  dpsidt(7) = 27*(xi(2)*xi(3)*dxidt(1) + xi(3)*xi(1)*dxidt(2) +...
                      xi(1)*xi(2)*dxidt(3)); 
      return
