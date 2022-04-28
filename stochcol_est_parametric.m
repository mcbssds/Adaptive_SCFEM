function perrests = stochcol_est_parametric(X_diff, errest2s, ...
    gridd_diff, list, rule_id)
%STOCHCOL_EST_PARAMETRIC indexwise parametric error estimator
%
% perrests = stochcol_est_parametric(X_diff, errest2s, ...
%    gridd_diff, list, rule_id)
%
%   TIFISS function: FX 28 November 2019
% Copyright (c) 2019 F. Xu

N = size(X_diff,1);
perrests = zeros(1, N);
for k = 1:size(X_diff,1)
    X_new = X_diff(k, :);
    perrests(k) = stochcol_interpost(X_new, errest2s, gridd_diff, ...
        list, rule_id);
end