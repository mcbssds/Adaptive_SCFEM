function KL_DATA = stoch_kl(input,expansion_type)
%STOCH_KL set up coefficients and their gradient of KL expansion
%  KL_DATA = stoch_kl(input,expansion_type);
%   SIFISS function: FX 8 March 2019
% Copyright (c) 2019 F. Xu

if expansion_type == 1
    norv2d = input(1);
    ax = input(2);
    ay = input(3);
    correl_x = input(4);
    correl_y = input(5);
    sigma = input(6);
    % set anonymous functions
    oddfun  = @(z,a,c) c-z*tan(a*z);
    evenfun = @(z,a,c) z+c*tan(a*z);
    % set the number of r.v. in each 1d direction
    norv1d = norv2d;
    % invert correlation lengths
    cx = 1/correl_x;
    cy = 1/correl_y;
    % initialisation
    omega_x = zeros(1, norv1d);
    omega_y = zeros(1, norv1d);
    lambda_x = zeros(1, norv1d);
    lambda_y = zeros(1, norv1d);
    alpha_x = zeros(1, norv1d);
    alpha_y = zeros(1, norv1d);
    ef_x = cell(1, norv1d);
    ef_y = cell(1, norv1d);
    for n = 1:norv1d
        if n==1
            x0 = [0, pi/(2*ax) - 1.0e-08];
            y0 = [0, pi/(2*ay) - 1.0e-08];
        else
            k = floor(n/2);
            x0 = [(2*k-1)*pi/(2*ax) + 1.0e-08, (2*k+1)*pi/(2*ax) - 1.0e-08];
            y0 = [(2*k-1)*pi/(2*ay) + 1.0e-08, (2*k+1)*pi/(2*ay) - 1.0e-08];
        end
        if mod(n,2) == 1
            omega_x(n) = fzero(@(z,a,c) oddfun(z,ax,cx), x0);
            omega_y(n) = fzero(@(z,a,c) oddfun(z,ay,cy), y0);
        else
            omega_x(n) = fzero(@(z,a,c) evenfun(z,ax,cx), x0);
            omega_y(n) = fzero(@(z,a,c) evenfun(z,ay,cy), y0);
        end
        lambda_x(n) = 2*cx/(omega_x(n)^2 + cx^2);
        lambda_y(n) = 2*cy/(omega_y(n)^2 + cy^2);
        alpha_x(n) = 1/sqrt(ax + (-1)^(n-1)*sin(2*ax*omega_x(n))/(2*omega_x(n)));
        alpha_y(n) = 1/sqrt(ay + (-1)^(n-1)*sin(2*ay*omega_y(n))/(2*omega_y(n)));
        if mod(n,2) == 1
            ef_x{n} = @(x) alpha_x(n)*(cos(omega_x(n)*x));
            ef_x_dx{n} = @(x) -alpha_x(n)*omega_x(n)*(sin(omega_x(n)*x));
            ef_y{n} = @(y) alpha_y(n)*(cos(omega_y(n)*y));
            ef_y_dy{n} = @(y) -alpha_y(n)*omega_y(n)*(sin(omega_y(n)*y));
        else
            ef_x{n} = @(x) alpha_x(n)*(sin(omega_x(n)*x));
            ef_x_dx{n} = @(x) alpha_x(n)*omega_x(n)*(cos(omega_x(n)*x));
            ef_y{n} = @(y) alpha_y(n)*(sin(omega_y(n)*y));
            ef_y_dy{n} = @(y) alpha_y(n)*omega_y(n)*(cos(omega_y(n)*y));
        end
    end
    % find 2d eigenvalues, sort them out keeping track of their 'directional' indices
    eigaux_2d = lambda_x' * lambda_y;
    [eigaux_2d_sorted, ind_sorted] = sort(eigaux_2d(:)','descend');
    [indaux_x, indaux_y] = ind2sub([norv1d norv1d],ind_sorted);
    % assign values to the fields of the structure KL_DATA
    KL_DATA.coeff{1} = @(x, y) ones(size(x)); % a_0
    KL_DATA.gradcoeff{1,1}= @(x, y) zeros(size(x)); % a_0_dx
    KL_DATA.gradcoeff{1,2}= @(x, y) zeros(size(y)); % a_0_dy
    for m = 1:norv2d
        ef_x_temp = ef_x{indaux_x(m)};
        ef_y_temp = ef_y{indaux_y(m)};
        ef_x_dx_temp = ef_x_dx{indaux_x(m)};
        ef_y_dy_temp = ef_y_dy{indaux_y(m)};
        ev2d_temp = eigaux_2d_sorted(m);
        KL_DATA.coeff{m+1,1} = @(x, y) sigma*sqrt(ev2d_temp)*ef_x_temp(x).*ef_y_temp(y); % a_m
        KL_DATA.gradcoeff{m+1,1} = @(x, y) sigma*sqrt(ev2d_temp)*ef_x_dx_temp(x).*ef_y_temp(y); % a_m_dx
        KL_DATA.gradcoeff{m+1,2} = @(x, y) sigma*sqrt(ev2d_temp)*ef_x_temp(x).*ef_y_dy_temp(y); % a_m_dy
    end
elseif expansion_type == 2 % Eigel
    norv = input(1);
    alpha_bar = input(2);
    sigma_tilde = input(3);
    KL_DATA.coeff{1} = @(x, y) ones(size(x)); % a_0
    KL_DATA.gradcoeff{1,1}= @(x, y) zeros(size(x)); % a_0_dx
    KL_DATA.gradcoeff{1,2}= @(x, y) zeros(size(y)); % a_0_dy
    for m = 1:norv
        km = floor(-0.5e0+sqrt(0.25e0+2*m));
        beta_x = m -km*(km + 1)/2;
        beta_y = km -beta_x;
        KL_DATA.coeff{m+1,1} = @(x, y) alpha_bar/(m^(sigma_tilde))*cos(2*pi*beta_x*x).*cos(2*pi*beta_y*y);
        KL_DATA.gradcoeff{m+1,1} = @(x, y) -2*pi*beta_x*alpha_bar/(m^(sigma_tilde))*sin(2*pi*beta_x*x).*cos(2*pi*beta_y*y);
        KL_DATA.gradcoeff{m+1,2} = @(x, y) -2*pi*beta_y*alpha_bar/(m^(sigma_tilde))*cos(2*pi*beta_x*x).*sin(2*pi*beta_y*y);
    end
end
end