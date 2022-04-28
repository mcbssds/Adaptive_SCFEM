function [x, y, xy, evt, bound, eboundt] = ell_domain_unstructured_fa
%ELL_DOMAIN_UNSTRUCTURED_FA  L-shape unstructured grid generator
%  ell_domain_unstructured_fa
%  input from keyboard by the default function:
%     h                        the mesh size
%     mesh_density_ratio       the relative mesh density close the corner,
%                              smaller mesh_density_ratio leads to finer
%                              grids close corner
%
% This is a copy of the original TIFISS function EEL_DOMAIN_UNSTRUCTURED 
% (DJS; 12 February 2016). The only difference here is in 
% - the default mesh size 'h'
% - the default 'mesh_density_ratio' 
% which are good choices for a coarse generation of the mesh.
%
% Function(s) called: distmesh2d (DISTMESH generation)
%                     find_boundary_elements
%                     mesh_regularity_check
%
% Grid defining data is saved to the file: datafiles/ell_grid.mat
%
%    TIFISS function: LR; 22 June 2018
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

global mesh_density_ratio   % called in the function ell_fh
% set DISTMESH parameters
h=default('Target mesh size (default 0.2)',0.2);
mesh_density_ratio = default('mesh density near corner (default 1)',1);
%it_number = 300;  % upper limit on number of mesh iterations
fprintf('Generating triangulation...');
[xy,evt]= distmesh2d(@ell_fd, @ell_fh, h, [-1,-1;1,1],[]);
fprintf('done\n');

%% get boundary nodes and edges
[eboundt,bound]=find_boundary_elements(evt);
%% check the mesh regularity
[xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);

%% reset boundary nodes
neboundt=length(eboundt(:,1));
for i=1:neboundt
    iel=eboundt(i,1);
    nodes=evt(iel,:);
    if eboundt(i,2)==1, nodeb=evt(iel,[2,3]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==2, nodeb=evt(iel,[3,1]); xybd=xy(nodeb,:); end
    if eboundt(i,2)==3, nodeb=evt(iel,[1,2]); xybd=xy(nodeb,:); end
    h=norm(xybd(1,:)-xybd(2,:),2);
    % vertical bound or horizontal bound
    if abs(xybd(1,2)-xybd(2,2))>abs(xybd(1,1)-xybd(2,1))
       if xybd(1,1)<-0.5, xy(nodeb,1)=-1;
       elseif xybd(1,1)<0.5, xy(nodeb,1)=0;
       else xy(nodeb,1)=1;
       end
    else
       if xybd(1,2)<-0.5, xy(nodeb,2)=-1;
       elseif xybd(1,2)<0.5, xy(nodeb,2)=0;
       else  xy(nodeb,2)=1;
       end
    end        
end

x = xy(:,1); y = xy(:,2);

gohome; cd datafiles
save ell_grid.mat xy evt bound eboundt x y

end % end function