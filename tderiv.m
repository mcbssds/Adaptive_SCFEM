function [jac,invjac,phi,dphidx,dphidy] = tderiv(s,t,xl,yl)
%TDERIV  evaluates derivatives of linear shape functions
%   [jac,invjac,phi,dphidx,dphidy] = tderiv(s,t,xl,yl);
%   input
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%   output
%          jac       elementwise jacobian (evaluated at (s,t))
%          invjac    elementwise inverse of jacobian
%          phi       elementwise shape functions
%          dphidx    x derivatives of phi
%          dphidy    y derivatives of phi
%
%    PIFISS function: DJS; 31 January 2007. 
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
      nel=length(xl(:,1));
      zero_v = zeros(nel,1);
      one_v = ones(nel,1);
%
% evaluate shape functions 
     [xi,dxids,dxidt] = tshape(s,t);
%
      dxds = zero_v;
      dxdt = zero_v;
      dyds = zero_v;
      dydt = zero_v;
	   jac = zero_v;
    invjac = zero_v; 
%
      for ivtx = 1:3
         dxds(:) = dxds(:) + xl(:,ivtx) .* one_v*dxids(ivtx);
         dxdt(:) = dxdt(:) + xl(:,ivtx) .* one_v*dxidt(ivtx);
         dyds(:) = dyds(:) + yl(:,ivtx) .* one_v*dxids(ivtx);
         dydt(:) = dydt(:) + yl(:,ivtx) .* one_v*dxidt(ivtx);
      end
%
      jac(:) = dxds(:).*dydt(:) - dxdt(:).*dyds(:);
      invjac(:) = one_v ./ jac(:);
%
      for ivtx = 1:3
         phi(:,ivtx) = xi(ivtx)*one_v;
		 dphidx(:,ivtx) = ( dxids(:,ivtx).*dydt(:) - dxidt(:,ivtx).*dyds(:));
         dphidy(:,ivtx) = (-dxids(:,ivtx).*dxdt(:) + dxidt(:,ivtx).*dxds(:));
      end
      return