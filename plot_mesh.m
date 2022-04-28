function plot_mesh(evt,xy,plot_title,iplot_elem,iplot_vtx,int,bound)
%PLOT_MESH plots mesh details using colors
%
% plot_mesh(evt,xy,plot_title,iplot_elem,iplot_vtx,int,bound)
%
% input:
%            evt      element mapping matrix
%             xy      vertex coordinate vector 
%     plot_title      (optional) figure's title
%     iplot_elem      (optional) element's number plotting switch
%      iplot_vtx      (optional) vertices's number plotting switch
%            int      (optional) vector of interior nodes  
%          bound      (optional) vector of boundary nodes
%
% The function plots the input mesh (given by evt and xy).
% It is also possible to display a given title, individual elements' numbers, 
% vertices' numbers, and blue squares ans red circles for each interanl and 
% boundary node, respectively.
%
% Default mesh colour is blue
%
% ------------------------------------------------------------------------
% Example 1. Calling without any optional input:
%    plot_mesh(evt,xy);
%
% Example 2. Calling with title, elements' numbers, and internal nodes:
%    plot_mesh(evt,xy,plot_title,1,[],int);
%
% Example 3. Calling with no title, with vertices' numbers, and both
% internal and boundary nodes:
%    plot_mesh(evt,xy,[],[],1,int,bound);
%
% Example 4:
%    xy = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5; 0.0 0.0];
%    evt = [1 5 4; 2 5 1; 3 5 2; 4 5 3];
%    int = [5];
%    bound = [1 2 3 4];
%    plot_mesh(evt,xy,'Mesh example',1,1,int,bound);
% ------------------------------------------------------------------------
%
% See also PLOTCOLOR_MARKEDELEM
%
%    TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi
  
  iplot_intbound = 1;
  if nargin < 7
     bound = [];                  % boundary nodes empty vector
     if nargin < 6
        int = [];                 % interior nodes empty vector
        iplot_intbound = 0;       % neither squares nor circles around vertices 
        if nargin < 5
           iplot_vtx = 0;         % no vertices' numbers
           if nargin < 4
              iplot_elem = 0;     % no elements' numbers
              if nargin < 3
                 plot_title = ''; % no title
                 if nargin < 2
                    error('Insufficient input number!');
                 end
              end
           end 
        end
     end
  end
            
  nel  = size(evt,1);   % number of elements
  nvtx = size(xy,1);    % number of nodes  
  rnxy = xy(int,:);     % internal nodes
  bnxy = xy(bound,:);   % boundary nodes

% -------------------------------------------------------------      
% Adjacency matrix
% -------------------------------------------------------------    
% This is a sparse adjacency matrix for the mesh given by evt matrix
  adj = sparse(nvtx,nvtx);
  evtt = [evt,evt(:,1)];
  for j = 1:3
      ncol1 = evtt(:,j);
      ncol2 = evtt(:,j+1);
      adj = adj + sparse(ncol1,ncol2,1,nvtx,nvtx);
  end

% -------------------------------------------------------------    
% Plot
% -------------------------------------------------------------    
  figure;
  hold on;
  gplot(adj,xy,'b');
 
% Display element's numbers  
  if iplot_elem
      xl_v = zeros(nel,3); 
      yl_v = zeros(nel,3); 
      for ivtx = 1:3
          xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all nodes
          yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all nodes
      end
      % Element's centroid coordinates
      xyc(:,1) = sum(xl_v,2) / 3;
      xyc(:,2) = sum(yl_v,2) / 3;
      %
      % Write element's number
      elem = (1:nel)';
      elenum = int2str(elem);      
      text(xyc(:,1),xyc(:,2),elenum,'Color','blue','Fontsize',12);
  end
  
% Display vertices' numbers
  if iplot_vtx      
      vertices = (1:nvtx)';
      vtxnum = int2str(vertices);
      text(xy(:,1),xy(:,2),vtxnum,'Color','black','Fontsize',14);
  end   
  
% Plotting a blue square and a red circle on each interior and boundary node, respectively
  if iplot_intbound
      plot(rnxy(:,1),rnxy(:,2),'bs');
      plot(bnxy(:,1),bnxy(:,2),'ro'); 
  end  
  
  axis('square');  
  title(plot_title,'Fontsize',17);
  set(gca,'FontSize',17);
  
% Axis ticks
  minx = min(xy(:,1)); maxx = max(xy(:,1)); midpx = 0.5*(minx+maxx); midplx = 0.5*(minx+midpx); midprx = 0.5*(midpx+maxx);
  miny = min(xy(:,2)); maxy = max(xy(:,2)); midpy = 0.5*(miny+maxy); midply = 0.5*(miny+midpy); midpry = 0.5*(midpy+maxy);
  xticks = [minx midplx midpx midprx maxx];
  yticks = [miny midply midpy midpry maxy];
  set(gca,'XTick',xticks,'YTick',yticks);
  hold off;

end  % end function