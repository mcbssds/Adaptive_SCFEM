function [p2xy,p2evt,p2bound]=p2grid(xy,evt,bound,eboundt,iplot)
%P2GRID   quadratic element grid generator
%[p2xy,p2evt,p2bound]=p2grid(xy,evt,bound,eboundt,iplot);
%   input
%          xy         vertex coordinate vector  
%          evt        element vertex mapping matrix
%          bound      boundary vertex vector
%          eboundt    element boundary mapping matrix
%          iplot      grid plotting switch
%   output
%          p2xy       P2 nodal coordinate vector
%          p2evt      P2 element node mapping matrix
%          p2bound    P2 boundary node vector
%    TIFISS function: DJS; 12 January 2016.
% Copyright (c) 2016 D.J. Silvester, Q. Liao
if nargin<5, iplot=1; end
nel=length(evt(:,1));
tedge=[];
for i=1:nel
   tedge=[tedge;evt(i,1),evt(i,2);evt(i,2),evt(i,3);evt(i,3),evt(i,1)];
end

%fprintf('start define edge')
ne=length(tedge(:,1));
for i=1:ne
    t=find((tedge(i,1)==tedge(1:(i-1),1)|(tedge(i,1)==tedge(1:(i-1),2)))&((tedge(i,2)==tedge(1:(i-1),1))|(tedge(i,2)==tedge(1:(i-1),2))) );
    if length(t)>0
       tedge(i,1)=nan;
    end
end
%fprintf('edges are defined\n')
%fprintf('start finding\n')
t=find(tedge(:,1)>0);
tedge=tedge(t,:);
ne=length(tedge);
eve=zeros(nel,3);
for i=1:nel
    eve(i,1)=find((evt(i,2)==tedge(:,1)|evt(i,2)==tedge(:,2))&(evt(i,3)==tedge(:,1)|evt(i,3)==tedge(:,2)));
    eve(i,2)=find((evt(i,1)==tedge(:,1)|evt(i,1)==tedge(:,2))&(evt(i,3)==tedge(:,1)|evt(i,3)==tedge(:,2)));
    eve(i,3)=find((evt(i,1)==tedge(:,1)|evt(i,1)==tedge(:,2))&(evt(i,2)==tedge(:,1)|evt(i,2)==tedge(:,2)));
end

middle_xy=zeros(ne,2);
n_node=length(xy(:,1));
for i=1:ne
    middle_xy(i,1)=(xy(tedge(i,1),1)+xy(tedge(i,2),1))/2;
    middle_xy(i,2)=(xy(tedge(i,1),2)+xy(tedge(i,2),2))/2;
end
p2xy=[xy;middle_xy];
p2evt=[evt,eve+n_node];

n_eboundt=length(eboundt);
p2bound=bound;
for i=1:n_eboundt
    p2bound=[p2bound;p2evt(eboundt(i,1),3+eboundt(i,2))];
end
nbound=length(p2bound);
for i=1:nbound-1
    for j=i+1:nbound
        if p2bound(j)==p2bound(i)
            p2bound(j)=nan;
        end
    end
end
t=find(p2bound>0);
p2bound=p2bound(t);
        
if (iplot)
% close the P1 grid figure
if ishandle(10), close 10, end
% plot the generated P2 grid
    nvtx=length(xy(:,1));
    adj=sparse(nvtx,nvtx);
    for i=1:nel
	adj(p2evt(i,1),p2evt(i,2)) =1;
  	adj(p2evt(i,2),p2evt(i,3)) =1;
	adj(p2evt(i,3),p2evt(i,1)) =1;
    end
    figure(11)
    gplot(adj,p2xy,'b')
    axis('equal')
    hold on
    plot(p2xy(:,1),p2xy(:,2),'ro')
    axis('off')
    xybd=p2xy(p2bound,:);
    plot(xybd(:,1),xybd(:,2),'ko')
    hold off
    title('P2 finite element subdivision');
    drawnow
end
return
