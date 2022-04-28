function sc_plotsol_p1(dom_type,sol,var,evt,xy)
%SC_PLOTSOL_P1 plots SC solution and variance for P1 approximation
%
% sc_plotsol_p1(dom_type,sol,var,evt,xy)
%
% input: 
%        dom_type      domain type
%             sol      nodal stochatic FE solution vector
%             var      nodal variance of the solution
%             evt      element mapping matrix
%              xy      vertex coordinate vector  
%
% NOTE that the mean-field and the variance are interpolated on a square grid 
% [X,Y] in order to plot isolines (contour plot); Matlab does not provide a 
% countour function for mesh-based functions.
%
%   TIFISS function: LR; 22 June 2018; DJS, 28 January 2021
% Copyright (c) 2018 A. Bespalov, L. Rocchi
    
  nvtx = size(xy,1);    % Number of vertices
  nel  = size(evt,1);   % Number of elements
 
% Refine grid and get cartesian product mesh
  npoints = 150; 
  x = linspace(min(xy(:,1)),max(xy(:,1)),npoints);
  y = linspace(min(xy(:,2)),max(xy(:,2)),npoints); 
  [X,Y] = meshgrid(x,y);
  
% -----------------------------------------------------------------------------  
% Expectation and variance of the solution
% -----------------------------------------------------------------------------
  expsol = griddata(xy(:,1),xy(:,2),sol(1:nvtx),X,Y);
  varsol = griddata(xy(:,1),xy(:,2),var(1:nvtx),X,Y);
% Fix to zero (eventual very small) negative values of the variance
% (due to griddata interpolation)
  varsol(varsol<0.0) = 0.0;
  
% Get rid of points outside the domain for X-Y grid
  [expsol,varsol] = recover_domain(dom_type,expsol,varsol,X,Y);  
 
% Plot mean value and variance
  title1 = 'Expectation';
  title2 = 'Variance';
  plot_stuff(X,Y,expsol,varsol,sol(1:nvtx),var,xy,evt,dom_type,title1,title2);
 

end % end function  


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function [func1,func2] = recover_domain(dom_type,func1,func2,X,Y)
% Get rid of those points of the cartesian grid X-Y outside the domain
  
% Nothing to do for convex domains (square):
% - dom_type == 1, (0,1)^2
  if dom_type == 2
      % L-shaped domain
      func1(X<0 & Y<0) = nan;
      func2(X<0 & Y<0) = nan;
  end
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function outline_domain(dom_type)
% Outline the corresponding domain 
  if dom_type == 1
      unitsquare;
  elseif dom_type == 2
      ellx;
  elseif dom_type==11
     bespokesquare
  end
  
end % end child function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function plot_stuff(X,Y,func1,func2,expmesh,varmesh,xy,evt,dom_type,title1,title2)
% In this plot function both mean value and variance are plotted using trimesh 

  fontSize = 12;

% Plot functions 
  figure%(fig);
  
  subplot(221)
  contour(X,Y,func1,20);
  axis square;  axis off; 
  outline_domain(dom_type);
  title(title1);
  
  subplot(222)
  trimesh(evt,xy(:,1),xy(:,2),expmesh); %mesh(X,Y,func1);          
  axis square;
  view(330,30);
  
  subplot(223)
  contour(X,Y,func2,20);    
  axis square;  axis off; 
  outline_domain(dom_type);
  title(title2);
   
  subplot(224)
  trimesh(evt,xy(:,1),xy(:,2),varmesh); % mesh(X,Y,func2);          
  axis square;
  view(330,30);
    
  set(findall(gcf,'-property','Fontsize'),'Fontsize',fontSize);
  
end % end child function


