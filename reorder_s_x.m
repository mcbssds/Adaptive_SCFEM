function [s,t] = reorder_s_x(sigpt,evt,edge_number)
%REORDER_S_X legacy code
% [s,t] = reorder_s_x(sigpt,evt,edge_number);
%   input
%          sigpt           1D Gauss Point
%          evt             element mapping matrix
%          edge_number     edge number
%   output
%          s,t             location of Gauss point on reference element
%   TIFISS function: QL; 14 April 2011.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

nel=length(evt(:,1));
s=zeros(nel,1);
t=zeros(nel,1);
if edge_number==1
    for i=1:nel
        if evt(i,2)<evt(i,3)
            s(i)=(1-sigpt)/2;
            t(i)=1-s(i);
        else
            s(i)=(1+sigpt)/2;
            t(i)=1-s(i);
        end
    end
elseif edge_number==2
    for i=1:nel
        if evt(i,3)<evt(i,1)
            s(i)=0;
            t(i)=(1-sigpt)/2;
        else
            s(i)=0;
            t(i)=(1+sigpt)/2;
        end
    end
elseif edge_number==3
    for i=1:nel
        if evt(i,1)<evt(i,2)
            s(i)=(1+sigpt)/2;
            t(i)=0;
        else
            s(i)=(1-sigpt)/2;
            t(i)=0;
        end
    end
end
            
            
