function [aa, a] = stochcol_diffusion_coeff_fun(KL_DATA, rf_type, M)

% function handle for KL expansion
a = @(x1, x2, yy) KL_DATA.coeff{1}(x1, x2);
for m = 1:M
    a = @(x1, x2, yy) a(x1, x2, yy) + KL_DATA.coeff{m+1}(x1, x2)*yy(m);
end
% function handles for exp(a), a^2 and a
if rf_type == 1
    aa = @(x1, x2, yy) exp(a(x1, x2, yy));
elseif rf_type == 2
    aa = @(x1, x2, yy) a(x1, x2, yy).^2;
elseif rf_type == 3
    aa = @(x1, x2, yy) a(x1, x2, yy);
end