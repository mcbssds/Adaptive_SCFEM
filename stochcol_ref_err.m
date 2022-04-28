function ref_err = stochcol_ref_err(sols, paras_sg, paras_fem, ...
    sols_ref, paras_sg_ref, paras_fem_ref, polys, G_ref, A_mean)
%STOCHCOL_REF_ERR compute reference error
%
% ref_err = stochcol_ref_err(sols, paras_sg, paras_fem, ...
%    sols_ref, paras_sg_ref, paras_fem_ref, polys, G_ref, A_mean)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

gridd          = paras_sg{4};
clincombiset   = paras_sg{5};
indlincombiset = paras_sg{6};
coords_ref = paras_sg_ref{end};
[N_ref, ~] = size(coords_ref);
xy        = paras_fem{1};
xy_ref    = paras_fem_ref{1};
% interpolate the approximation in reference collocation points and
% reference finite element grid
sols_interp = zeros(size(sols_ref));
for k = 1:N_ref
    LL = stochcol_getinterpolant_2(gridd, ...
        clincombiset, indlincombiset, coords_ref(k,:), polys);
    uyy = sols*LL;
    uyy_interp = griddata(xy(:,1), xy(:,2), uyy, xy_ref(:,1), xy_ref(:,2));
    sols_interp(:,k) = uyy_interp;
end
sols_diff = sols_ref - sols_interp;
ref_err = sqrt(sum(dot(G_ref, sols_diff' * A_mean * sols_diff)));

