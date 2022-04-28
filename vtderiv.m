function [jac,invjac,phi,dphidx,dphidy] = vtderiv(s,t,xl,yl)
%VTDERIV  derivatives of linear shape functions for vector s, t
%   [jac,invjac,phi,dphidx,dphidy] = vtderiv(s,t,xl,yl);
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
%    TIFISS function: QL; 17 April 2011.
% Copyright (c) 2011 D.J. Silvester and Q. Liao
      nel=length(xl(:,1));
      zero_v = zeros(nel,1);
      one_v = ones(nel,1);
%
% evaluate shape functions 
     [xi,dxids,dxidt] = vtshape(s,t);
%
      dxds = zero_v;
      dxdt = zero_v;
      dyds = zero_v;
      dydt = zero_v;
	   jac = zero_v;
    invjac = zero_v; 
%
      for ivtx = 1:3
         dxds(:) = dxds(:) + xl(:,ivtx) .* one_v.*dxids(:,ivtx);
         dxdt(:) = dxdt(:) + xl(:,ivtx) .* one_v.*dxidt(:,ivtx);
         dyds(:) = dyds(:) + yl(:,ivtx) .* one_v.*dxids(:,ivtx);
         dydt(:) = dydt(:) + yl(:,ivtx) .* one_v.*dxidt(:,ivtx);
      end
%
      jac(:) = dxds(:).*dydt(:) - dxdt(:).*dyds(:);
      invjac(:) = one_v ./ jac(:);
%
      for ivtx = 1:3
         phi(:,ivtx) = xi(:,ivtx).*one_v;
		 dphidx(:,ivtx) = ( dxids(:,ivtx).*dydt(:) - dxidt(:,ivtx).*dyds(:));
         dphidy(:,ivtx) = (-dxids(:,ivtx).*dxdt(:) + dxidt(:,ivtx).*dxds(:));
      end
      return