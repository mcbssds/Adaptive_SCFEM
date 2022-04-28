function [elerr_p, fe, ae] = diffpost_p2_with_p4(xy, evt, eex, tve, ...
    els, eboundt, p2sol, coeff_fun, rhs_fun)
%DIFFPOST_P2_WITH_P4 local Poisson error estimator for P2 with P4 correction
%
%   [elerr_p, fe, ae] = diffpost_p2_with_p4(xy, evt, eex, tve, ...
%    els, eboundt, p2sol, coeff_fun, rhs_fun)
%
%   input
%          xy         vertex coordinate vector
%          evt        element mapping matrix
%          eex        presorted element connectivity array
%          tve        presorted edge location array
%          els        element edge size array
%          eboundt    element edge boundary matrix
%          p2sol      P2 solution for diffusion problem
%        coeff_fun    function handle for the diffusion coefficient
%          rhs_fun    function handle for the right-hand side
%   output
%          elerr_p    elementwise error estimate
%          fe         elementwise rhs vectors
%          ae         elementwise Poisson problem matrices
%
% calls function:
%            tgauss
%            triangular_gausspoints
%            gausspoints_oned
%            tderiv
%            vtqarderiv
%            vtqderiv_2
%            p2fluxjmps
%            vtqarshape
%
%   TIFISS function: DJS; 15 January 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao

fprintf('fast computation of P2 local error estimator...\n')
x = xy(:, 1); 
y = xy(:, 2);
nel = length(evt(:, 1));
elerr_p = zeros(nel, 1);

% define reduced P4 approximation
node_basis = [7; 8; 9; 10; 11; 12; 13; 14; 15];
nnode = length(node_basis);

%% set up 19 point Gauss Rule
nngpt = 19; 
[s, t, wt] = triangular_gausspoints(nngpt); %% nngpt=73

% inner loop over elements
for ivtx = 1:3
    xl_v(:,ivtx) = x(evt(:,ivtx));
    yl_v(:,ivtx) = y(evt(:,ivtx));
end
for ivtx=1:6
    sl_v(:,ivtx) = p2sol(evt(:,ivtx));
end
ae = zeros(nel, nnode, nnode); 
elerr = zeros(nnode, nel);
fe = zeros(nel, nnode);
% loop over Gauss points
for igpt = 1:nngpt
    sigpt = s(igpt); 
    tigpt = t(igpt); 
    wght  = wt(igpt);
    % evaluate derivatives etc
    [jac_v, invjac_v, phi_v, dphidx_v, dphidy_v] = tderiv(sigpt, tigpt, ...
        xl_v, yl_v);
    [qar_v, dqardx_v, dqardy_v] = vtqarderiv(sigpt, tigpt, xl_v, yl_v);
    rhs_v = tgauss(sigpt, tigpt, xl_v, yl_v, rhs_fun);
    coeff = tgauss(sigpt, tigpt, xl_v, yl_v, coeff_fun);
    [dpsi2dx2, dpsi2dy2] = vtqderiv_2(sigpt.*ones(nel, 1), ...
        tigpt.*ones(nel,1), xl_v, yl_v);
    duh2dx2 = zeros(nel, 1);
    duh2dy2 = zeros(nel, 1);
    for i = 1:6
        duh2dx2(:) = duh2dx2(:) + sl_v(:, i).*dpsi2dx2(:, i).*invjac_v.^2;
        duh2dy2(:) = duh2dy2(:) + sl_v(:, i).*dpsi2dy2(:, i).*invjac_v.^2;
    end
    for j = 1:nnode
        jj = node_basis(j);
        for i = 1:nnode
            ii = node_basis(i);
            ae(:, i, j) = ae(:,i,j) + ...
                wght*dqardx_v(:, ii).*dqardx_v(:, jj).*invjac_v(:);
            ae(:, i, j) = ae(:,i,j) + ...
                wght*dqardy_v(:, ii).*dqardy_v(:, jj).*invjac_v(:);
        end
        fe(:, j) = fe(:, j) + wght*(coeff(:).*duh2dx2 + ...
            coeff(:).*duh2dy2 + rhs_v).*qar_v(:, jj).*jac_v(:);
    end
    % end of Gauss point loop
end
%
% compute flux jumps
% set up one-dimensional Gauss Rule
nngpt = 3; 
[oneg,onew] = gausspoints_oned(nngpt); %%% nngpt=7
for ig=1:nngpt
    sigpt = oneg(ig); 
    wght  = onew(ig);
    jmp = p2fluxjmps(p2sol, eex, tve, xy, evt, eboundt, sigpt, coeff_fun);
    for ee=1:nnode
        eei=node_basis(ee);
        if eei==2| eei==7| eei==4| eei==8| eei==3
            [s,t] = reorder_s_x(sigpt,evt,1);
            [qar_v,dqardx_v,dqardy_v] = vtqarderiv(s,t,xl_v,yl_v);
            fe(:,ee) = fe(:,ee) - wght*1/2*jmp(:,1).*qar_v(:,eei).* els(:,1)./2;
        end
        if eei==3| eei==9| eei==5| eei==10| eei==1
            [s,t] = reorder_s_x(sigpt,evt,2);
            [qar_v,dqardx_v,dqardy_v] = vtqarderiv(s,t,xl_v,yl_v);
            fe(:,ee) = fe(:,ee) - wght*1/2*jmp(:,2).*qar_v(:,eei).* els(:,2)./2;
        end
        if eei==1| eei==11| eei==6| eei==12| eei==2
            [s,t] = reorder_s_x(sigpt,evt,3);
            [qar_v,dqardx_v,dqardy_v] = vtqarderiv(s,t,xl_v,yl_v);
            fe(:,ee) = fe(:,ee) - wght*1/2*jmp(:,3).*qar_v(:,eei).* els(:,3)./2;
        end
    end
end


% solve for local estimate (sequential code) with middle point constraint
[qar,dqards,dqardt] = vtqarshape(1/3,1/3);
center_constrain_node = qar(1,node_basis);

%         for ielem = 1:nel
%		    elerr(:,ielem) = squeeze(ae(ielem,1:nnode,1:nnode))\(fe(ielem,1:nnode)');
%         end

% vectorized code
% LDLT factorization
nn=nnode;
dd=zeros(nel,nn); rr=zeros(nel,nn);
for kk=1:nn-1
    for pp=1:kk-1;
        rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,kk,pp);
    end
    dd(1:nel,kk)=ae(1:nel,kk,kk);
    for pp=1:kk-1;
        dd(1:nel,kk)= dd(1:nel,kk)- ae(1:nel,kk,pp).*rr(1:nel,pp);
    end
    for ii=kk+1:nn
        for pp=1:kk-1;
            ae(1:nel,ii,kk)=ae(1:nel,ii,kk)-ae(1:nel,ii,pp).*rr(1:nel,pp);
        end
        ae(1:nel,ii,kk)=ae(1:nel,ii,kk)./dd(1:nel,kk);
    end
end
for pp=1:nn-1;
    rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,nn,pp);
end
dd(1:nel,nn)=ae(1:nel,nn,nn);
for pp=1:nn-1;
    dd(1:nel,nn)= dd(1:nel,nn)- ae(1:nel,nn,pp).*rr(1:nel,pp);
end
% overwrite diagonal entries
for kk=1:nn
    ae(1:nel,kk,kk)= dd(1:nel,kk);
end
%  forward-backward substitutions ...
xx = element_lusolve(ae,fe);
elerr=xx';

for ivtx=1:nnode
    elerr_p(:) = elerr_p(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
elerr_p=sqrt(elerr_p);
fprintf('estimated energy error is %10.4e \n',norm(elerr_p,2))
return
