function jmp = p2fluxjmps(p2sol, eex, tve, xy, evt, eboundt, sigpt, c_fun)
%P2FLUXJMPS  vectorised flux jumps for triangular P2 grid
%
%  jmp = p2fluxjmps(p2sol, eex, tve, xy, evt, eboundt, sigpt, c_fun)
%
%   input
%          p2sol      nodal solution vector
%            eex      presorted element connectivity array
%            tve      presorted edge location array
%             xy      vertex coordinate vector
%            evt      element mapping matrix
%        eboundt      element edge boundary matrix
%          sigpt      Gauss point in (-1,1)
%          c_fun      function handle for the diffusion coefficient
%
%   output
%            jmp      component elementwise edge flux jumps
%
% Function(s) called: reorder_s_x
%                     vtderiv
%                     vtqderiv
%                     tgauss
%
%   TIFISS function: DJS; 16 January 2016.
% Copyright (c) 2016 D.J. Silvester and Qifeng Liao
%
x = xy(:, 1);
y = xy(:, 2);
nvtx = length(x);
nel = length(evt(:, 1));
%fprintf('computing P2 flux jumps... ')
%
% initialise global matrices
jmp = zeros(nel, 3);
flux = zeros(nel, 4);
zero_v = zeros(nel, 1);
one_v = ones(nel, 1);
%
% inner loop over elements
for ivtx = 1:3
    xl_v(:, ivtx) = x(evt(:, ivtx));
    yl_v(:, ivtx) = y(evt(:, ivtx));
end
for ivtx = 1:6
    sl_v(:, ivtx) = p2sol(evt(:, ivtx));
end

% evaluate derivatives and normal flux on each edge in turn
% first edge
[s, t] = reorder_s_x(sigpt, evt, 1);
[jac, invjac, phi, dphidx, dphidy] = vtderiv(s, t, xl_v, yl_v);
[psi, dpsidx, dpsidy] = vtqderiv(s, t, xl_v, yl_v);
[coeff] = tgauss(s, t, xl_v, yl_v, c_fun);
hx_v = xl_v(:,3) - xl_v(:,2);
hy_v = yl_v(:,3) - yl_v(:,2);
he_v = sqrt(hx_v.*hx_v+hy_v.*hy_v);
sx_v = hx_v./he_v;
sy_v = hy_v./he_v;  % unit tangential components
nx_v = sy_v;
ny_v = -sx_v;              % unit normal components
fx_v = zero_v;
for  ivtx = 1:6
    fx_v = fx_v + (coeff(:).*dpsidx(:,ivtx).*nx_v + ...
        coeff(:).*dpsidy(:,ivtx).*ny_v ).*invjac(:).*sl_v(:,ivtx);
end
flux(:, 1) = fx_v;

%
%% second edge
[s,t] = reorder_s_x(sigpt, evt, 2);
[jac, invjac, phi, dphidx, dphidy] = vtderiv(s, t, xl_v, yl_v);
[psi, dpsidx, dpsidy] = vtqderiv(s, t, xl_v, yl_v);
[coeff] = tgauss(s, t, xl_v, yl_v, c_fun);
hx_v = xl_v(:,1) - xl_v(:,3);
hy_v = yl_v(:,1) - yl_v(:,3);
he_v = sqrt(hx_v.*hx_v+hy_v.*hy_v);
sx_v = hx_v./he_v;
sy_v = hy_v./he_v;  % unit tangential components
nx_v = sy_v;
ny_v = -sx_v;              % unit normal components
fx_v = zero_v;
for  ivtx=1:6
    fx_v = fx_v + (coeff(:).*dpsidx(:,ivtx).*nx_v + ...
        coeff(:).*dpsidy(:,ivtx).*ny_v ).*invjac(:).*sl_v(:,ivtx);
end
flux(:,2)=fx_v;
%
%% third edge
[s, t] = reorder_s_x(sigpt, evt, 3);
[jac, invjac, phi, dphidx, dphidy] = vtderiv(s, t, xl_v, yl_v);
[psi, dpsidx, dpsidy] = vtqderiv(s, t, xl_v, yl_v);
[coeff] = tgauss(s, t, xl_v, yl_v, c_fun);
hx_v = xl_v(:,2)-xl_v(:,1);
hy_v = yl_v(:,2)-yl_v(:,1);
he_v = sqrt(hx_v.*hx_v + hy_v.*hy_v);
sx_v = hx_v./he_v;
sy_v = hy_v./he_v;  % unit tangential components
nx_v = sy_v;
ny_v = -sx_v;              % unit normal components
fx_v = zero_v;
for  ivtx=1:6
    fx_v = fx_v + (coeff(:).*dpsidx(:,ivtx).*nx_v + ...
        coeff(:).*dpsidy(:,ivtx).*ny_v ).*invjac(:).*sl_v(:,ivtx);
end
flux(:,3) = fx_v;

% add zero column for boundary jumps
flux(:,4) = zero_v;

% replace zero indices in  array tve by 4s
tvx = tve;
tvx(find(tve==0)) = 4;

%
% evaluate flux jump on each edge in turn
% A(sub2ind(size(A),ii,jj)) pulls out the entries of flux indexed by ii and jj
% first edge
jmp(:,1) = flux(:,1) + flux(sub2ind([nel,4],eex(:,1),tvx(:,1)) );
% second edge
jmp(:,2) = flux(:,2) + flux(sub2ind([nel,4],eex(:,2),tvx(:,2)) );
% third edge
jmp(:,3) = flux(:,3) + flux(sub2ind([nel,4],eex(:,3),tvx(:,3)) );

%
% remove Dirichlet boundary edge contributions
nbde = length(eboundt(:,1));
for k=1:nbde
    jmp(eboundt(k,1),eboundt(k,2))=0;
end
%
%      fprintf('done\n')
return
