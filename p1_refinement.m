function [rxy,revt,rbound,reboundt]=p1_refinement(xy,evt,bound,eboundt,iplot)
%P1_REFINEMENT uniformly refine the triangular mesh
%   input
%          xy         vertex coordinate vector  
%          evt        element vertex mapping matrix
%          bound      boundary vertex vector
%          eboundt    element boundary mapping matrix
%          iplot      grid plotting switch
%   output
%          rxy        refined vertex coordinate vector  
%          revt       refined element vertex mapping matrix
%          rbound     refined boundary vertex vector
%          reboundt   refined element boundary mapping matrix
%    TIFISS function: DJS; 12 January 2016.
% Copyright (c) 2016 D.J. Silvester, Qifeng Liao
if nargin<5, iplot=1; end
[p2xy,p2evt,p2bound]=p2grid(xy,evt,bound,eboundt,0);

rxy=p2xy;
nel=length(evt(:,1));
for i=1:nel
    revt((i-1)*4+1,:)=p2evt(i,[1,6,5]);
    revt((i-1)*4+2,:)=p2evt(i,[6,2,4]);
    revt((i-1)*4+3,:)=p2evt(i,[5,4,3]);
    revt((i-1)*4+4,:)=p2evt(i,[4,5,6]);
end
rbound=p2bound;
neb=length(eboundt(:,1));
for i=1:neb
    rnb=eboundt(i,1);
    reb=eboundt(i,2);
    if reb==1
        t1=2;
        t2=3;
    elseif reb==2
        t1=3;
        t2=1;
    elseif reb==3
        t1=1;
        t2=2;
    end
    reboundt((i-1)*2+1,1)=(rnb-1)*4+t1;
    reboundt((i-1)*2+1,2)=reb;
    reboundt((i-1)*2+2,1)=(rnb-1)*4+t2;
    reboundt((i-1)*2+2,2)=reb;
end

if (iplot)
% plotting the generated grid
    nvtx=length(rxy(:,1));
    rnel=length(revt(:,1));
	adj=sparse(nvtx,nvtx);
    for i=1:rnel
	adj(revt(i,1),revt(i,2)) =1;
	adj(revt(i,2),revt(i,3)) =1;
	adj(revt(i,3),revt(i,1)) =1;
    end
    figure(10)
    gplot(adj,rxy,'b')
    axis('square')
    hold on
    plot(rxy(:,1),rxy(:,2),'ro')
    rxybd=rxy(rbound,:);
    plot(rxybd(:,1),rxybd(:,2),'ko')
    hold off
    title('P1 finite element subdivision');
    drawnow
end