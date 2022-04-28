function [s,t,wt]=assemble_gausspoints(sd,td,rd,wo)
%ASSEMBLE_GAUSSPOINTS generates Gauss points from triangular coordinates
% [s,t,wt] = assemble_gausspoints(sd,td,rd,wo);
%  input
%      sd td rd       triangular coordinates
%      wo             weights
%  output
%      s              x-coordinate of the Gauss points
%      t              y-coordinate of the Gauss points
%      w              weights
%    TIFISS function: QL; 17 April 2011.
%  Copyright (c) 2011 D.J. Silvester, Qifeng Liao

% triangular coodinates are defined via vertex ordering:
%
%   v_2:(0,1)          sd: the triangular coordinate facing v_1
%      |\               
%      | \             td: the triangular coordinate facing v_2
%      |  \ 
%      |___\           rd: the triangular coordinate facing v_3
% v_3:(0,0) v_1:(1,0)

% note the three possible combinations
% 1: sd==td==rd;
% 2: td==rd, while sd~=td;
% 3: sd td rd are different numbers.


n=length(wo);
j=1;
for i=1:n
    % if all the triangular coordinates are same
    if sd(i)==td(i) 
        if td(i)~=rd(i),error('incorrect Gauss point ordering'), end
        wt(j)=wo(i);
        s(j)=sd(i);
        t(j)=td(i);
        j=j+1;
    % if two triangular coordinates are same (make sure they are td and sd)
    elseif td(i)==rd(i)
           wt(j)=wo(i);
           s(j)=sd(i);
           t(j)=td(i);
           j=j+1;
           wt(j)=wo(i);
           s(j)=td(i);
           t(j)=sd(i);
           j=j+1;
           wt(j)=wo(i);
           s(j)=td(i);
           t(j)=rd(i);
           j=j+1;
    % if none of the triangular coordinates is same
    else
        wt(j)=wo(i);
        s(j)=sd(i);
        t(j)=td(i);
        j=j+1;
        wt(j)=wo(i);
        s(j)=sd(i);
        t(j)=rd(i);       
        j=j+1;
        
        wt(j)=wo(i);
        s(j)=td(i);
        t(j)=sd(i);
        j=j+1;
        wt(j)=wo(i);
        s(j)=td(i);
        t(j)=rd(i);       
        j=j+1;
        
        wt(j)=wo(i);
        s(j)=rd(i);
        t(j)=sd(i);
        j=j+1;
        wt(j)=wo(i);
        s(j)=rd(i);
        t(j)=td(i);       
        j=j+1;     
    end
end
return
        
         