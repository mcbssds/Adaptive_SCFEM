function [aax1, aax2] = stochcol_diffusion_gradcoeff_fun(KL_DATA, rf_type, M, aa, a)

% gradients of KL expansion
ax1 = @(x1, x2, yy) KL_DATA.gradcoeff{1, 1}(x1, x2);
ax2 = @(x1, x2, yy) KL_DATA.gradcoeff{1, 2}(x1, x2);
for m = 1:M
    ax1 = @(x1, x2, yy) ax1(x1, x2, yy) + ...
        KL_DATA.gradcoeff{m+1, 1}(x1, x2)*yy(m);
    ax2 = @(x1, x2, yy) ax2(x1, x2, yy) + ...
        KL_DATA.gradcoeff{m+1, 2}(x1, x2)*yy(m);
end
if rf_type == 1
    aax1 = @(x1, x2, yy) aa(x1, x2, yy).*ax1(x1, x2, yy);
    aax2 = @(x1, x2, yy) aa(x1, x2, yy).*ax2(x1, x2, yy);
elseif rf_type == 2
    aax1 = @(x1, x2, yy) 2*a(x1, x2, yy).*ax1(x1, x2, yy);
    aax2 = @(x1, x2, yy) 2*a(x1, x2, yy).*ax2(x1, x2, yy);
elseif rf_type == 3
    aax1 = @(x1, x2, yy) ax1(x1, x2, yy);
    aax2 = @(x1, x2, yy) ax2(x1, x2, yy);
end