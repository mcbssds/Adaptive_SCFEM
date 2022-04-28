function [evt,eboundt]=p1grid(xy,mv,bound,mbound,iplot);
%P1GRID   linear element grid generator
%   [evt,eboundt]=p1grid(xy,mv,bound,mbound,iplot);
%   input
%          xy         vertex coordinate vector  
%          mv         Q2 macroelement mapping matrix
%          bound      boundary vertex vector
%          mbound     macroelement boundary edge vector
%          iplot      grid plotting switch
%   output
%          evt         element mapping matrix
%          eboundt     element boundary mapping matrix
%
%    PIFISS function: DJS; 12 January 2016.
% Copyright (c) 2007 C.E. Powell, D.J. Silvester
if nargin<5, iplot=1; end

xx=xy(:,1); yy=xy(:,2); nvtx=length(xx);
adj=sparse(nvtx,nvtx);
mel=length(mv(:,1)); nel=8*mel;
evt=zeros(nel,3);
%
%% loop over macroelements
for k=1:mel
% first element
ke=8*k-7;
evt(ke,1)=mv(k,9);
evt(ke,2)=mv(k,8);
evt(ke,3)=mv(k,1);
% second element
ke=8*k-6;
evt(ke,1)=mv(k,1);
evt(ke,2)=mv(k,5);
evt(ke,3)=mv(k,9);
% third element
ke=8*k-5;
evt(ke,1)=mv(k,9);
evt(ke,2)=mv(k,5);
evt(ke,3)=mv(k,2);
% fourth element
ke=8*k-4;
evt(ke,1)=mv(k,2);
evt(ke,2)=mv(k,6);
evt(ke,3)=mv(k,9);
% fifth element
ke=8*k-3;
evt(ke,1)=mv(k,9);
evt(ke,2)=mv(k,6);
evt(ke,3)=mv(k,3);
% sixth element
ke=8*k-2;
evt(ke,1)=mv(k,3);
evt(ke,2)=mv(k,7);
evt(ke,3)=mv(k,9);
% seventh element
ke=8*k-1;
evt(ke,1)=mv(k,9);
evt(ke,2)=mv(k,7);
evt(ke,3)=mv(k,4);
% eighth element
ke=8*k;
evt(ke,1)=mv(k,4);
evt(ke,2)=mv(k,8);
evt(ke,3)=mv(k,9);
end
%
%% define element edges
ect=1;
% bottom boundary edges
k1=find(mbound(:,2)==1)';
for k=mbound(k1)
   eboundt(ect,1)=8*k-6; eboundt(ect+1,1)=8*k-5; 
   eboundt(ect,2)=3    ; eboundt(ect+1,2)=1;
   ect=ect+2;
end
% right boundary edges
k2=find(mbound(:,2)==2)';
for k=mbound(k2)
   eboundt(ect,1)=8*k-4; eboundt(ect+1,1)=8*k-3; 
   eboundt(ect,2)=3    ; eboundt(ect+1,2)=1;
   ect=ect+2;
end
% top boundary edges
k3=find(mbound(:,2)==3)';
for k=mbound(k3)
   eboundt(ect,1)=8*k-2; eboundt(ect+1,1)=8*k-1; 
   eboundt(ect,2)=3    ; eboundt(ect+1,2)=1;
   ect=ect+2;
end
% left boundary edges
k4=find(mbound(:,2)==4)';
for k=mbound(k4)
   eboundt(ect,1)=8*k; eboundt(ect+1,1)=8*k-7; 
   eboundt(ect,2)=3    ; eboundt(ect+1,2)=1;
   ect=ect+2;
end
%

if (iplot)
% plot the generated grid
	adj=sparse(nvtx,nvtx);
    for i=1:nel
	adj(evt(i,1),evt(i,2)) =1;
	adj(evt(i,2),evt(i,3)) =1;
	adj(evt(i,3),evt(i,1)) =1;
    end
    figure(10)
    gplot(adj,xy,'b')
    axis('equal')
    hold on
    plot(xy(:,1),xy(:,2),'ro')
    xybd=xy(bound,:);
    plot(xybd(:,1),xybd(:,2),'ko')
    hold off
    title('P1 finite element subdivision');
    drawnow
end
return
