function perrest = stochcol_interpost(X_new, errest2s, gridd_diff, ...
    list, rule_id)
%STOCHCOL_INTERPOST compute parametric error
%
% perrest = stochcol_interpost(X_new, errest2s, gridd_diff, ...
%    list, rule_id)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

paras_sg_new = stochcol_sg(X_new, rule_id);
L_two_norm_new = stochcol_multilag_Ltwonorm(paras_sg_new, list);
gridd_new = paras_sg_new{4};
[~,IA, IB] = intersect(gridd_diff, gridd_new, 'rows');
perrest = errest2s(IA)*L_two_norm_new(IB);
