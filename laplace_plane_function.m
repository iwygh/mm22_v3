function [p0,q0,i0,Omega0] = laplace_plane_function(a_sample,planetdir,planetfile)
    mu_sun = 132712440041.93938; % km^3/s^2
    mu_sun = mu_sun*1e9; % m^3/s^2
    au = 1.495978707e11; % m
    a_sample = a_sample * au; % m
    n_sample = sqrt(mu_sun/a_sample^3); % rad/s
    cd(planetdir)
    tbl = readtable(planetfile);
    GM_solar_masses_mercury = tbl.GM_solar_masses(2);
    GM_solar_masses_venus   = tbl.GM_solar_masses(3);
    GM_solar_masses_earth   = tbl.GM_solar_masses(4);
    GM_solar_masses_mars    = tbl.GM_solar_masses(5);
    GM_solar_masses_jupiter = tbl.GM_solar_masses(6);
    GM_solar_masses_saturn  = tbl.GM_solar_masses(7);
    GM_solar_masses_uranus  = tbl.GM_solar_masses(8);
    GM_solar_masses_neptune = tbl.GM_solar_masses(9);
    GM_solar_masses_sun = 1 + GM_solar_masses_mercury + GM_solar_masses_venus + ...
                              GM_solar_masses_earth   + GM_solar_masses_mars;
    GM_solar_masses_jupiter = GM_solar_masses_jupiter / GM_solar_masses_sun;
    GM_solar_masses_saturn  = GM_solar_masses_saturn  / GM_solar_masses_sun;
    GM_solar_masses_uranus  = GM_solar_masses_uranus  / GM_solar_masses_sun;
    GM_solar_masses_neptune = GM_solar_masses_neptune / GM_solar_masses_sun;
    GM_solar_masses_sun     = GM_solar_masses_sun     / GM_solar_masses_sun;
    mu_jupiter = GM_solar_masses_jupiter * mu_sun;
    mu_saturn  = GM_solar_masses_saturn  * mu_sun;
    mu_uranus  = GM_solar_masses_uranus  * mu_sun;
    mu_neptune = GM_solar_masses_neptune * mu_sun;
    mu_sun     = GM_solar_masses_sun     * mu_sun;
    mu_planets = [mu_jupiter;
                  mu_saturn;
                  mu_uranus;
                  mu_neptune];
    a_planets = tbl.semimajor_axis_au(6:9) * au; % 6:9 means take the outer planets only
    i_planets = tbl.inclination_degrees(6:9) * pi/180;
    Omega_planets = tbl.longitude_of_node_degrees(6:9) * pi/180;
%     mu_planets = [0.00095479191521124;
%                   0.000285885672722241;
%                   4.36624373583127E-05;
%                   5.15138377262867E-05] * mu_sun;
%     a_planets = [5.20309819923364;
%                  9.58132546622039;
%                  19.2502696807627;
%                  30.2983527657524] * au;
%     i_planets = [1.30358696992995;
%                  2.48633767637901;
%                  0.770433821452291;
%                  1.76887985910468]; % deg
%     Omega_planets = [100.517119794776;
%                      113.59571858058;
%                      74.0913281240404;
%                      131.741931958958]; % deg
%     i_planets = i_planets*pi/180; % rad
%     Omega_planets = Omega_planets*pi/180; % rad
    n_planets = sqrt(mu_sun./a_planets.^3); % rad/s
    N = length(a_planets); % == 8
    B = zeros(N,N);
    alpha = zeros(N,N);
    alphabar = zeros(N,N);
    % eq 7.13
    b32_1_integrand = @(psi,alpha) cos(psi)./(1-2*alpha*cos(psi)+alpha.^2).^(3/2);
    b32_1_fun = @(alpha) integral(@(psi) b32_1_integrand(psi,alpha),0,2*pi) / pi;
    % eq 7.128, 7.129
    for j = 1:N
        for k = 1:N
           aj = a_planets(j);
           ak = a_planets(k);
           if aj > ak
               alpha(j,k) = ak/aj;
               alphabar(j,k) = 1;
           else
               alpha(j,k) = aj/ak;
               alphabar(j,k) = aj/ak;
           end
        end
    end
    % eq 7.134, 7.135
    warning('off','MATLAB:integral:NonFiniteValue')
    for j = 1:N
        for k = 1:N
            B(j,k) = 1/4 * mu_planets(k)/(mu_sun+mu_planets(j)) * ...
                n_planets(j) * alpha(j,k) * alphabar(j,k) * b32_1_fun(alpha(j,k));
        end
    end
    for j = 1:N
        B(j,j) = 0;
        for k = 1:N
            if k ~= j
                B(j,j) = B(j,j) - n_planets(j) * 1/4 * mu_planets(k) / ...
                    (mu_sun+mu_planets(j)) * alpha(j,k) * alphabar(j,k) * ...
                    b32_1_fun(alpha(j,k));
            end
        end
    end
    [VB,DB] = eig(B);
    f = diag(DB); % pg 301 below eq 7.138
    Ibar = VB; % pg 301 below eq 7.138
    p_planets = i_planets.*sin(Omega_planets); % eq 7.19
    q_planets = i_planets.*cos(Omega_planets); % eq 7.19
    T_singamma = Ibar\p_planets; % eq 7.47
    T_cosgamma = Ibar\q_planets; % eq 7.47
    T = sqrt(T_singamma.^2+T_cosgamma.^2);
    singamma = T_singamma./T;
    cosgamma = T_cosgamma./T;
    I = zeros(N,N);
    for i = 1:N
        for j = 1:N
            I(j,i) = Ibar(j,i)*T(i); % eq 7.41
        end
    end
    gamma = atan2(singamma,cosgamma);
    gamma = gamma + (gamma<0)*2*pi;
    B_vec = zeros(N,1);
    alpha_vec = zeros(N,1);
    alphabar_vec = zeros(N,1);
    % eq 7.128, 7.129
    for i = 1:N
        if a_planets(i) < a_sample
            alpha_vec(i) = a_planets(i)/a_sample;
            alphabar_vec(i) = 1;
        else
            alpha_vec(i) = a_sample/a_planets(i);
            alphabar_vec(i) = a_sample/a_planets(i);
        end
    end
    % eq 7.144
    for i = 1:N
        B_vec(i) = n_sample/4*mu_planets(i)/mu_sun*alpha_vec(i)*alphabar_vec(i) * ...
            b32_1_fun(alpha_vec(i));
    end
    % eq 7.143
    B_scalar = -sum(B_vec);
    % eq 7.76
    mu_vec = zeros(N,1); % not mu as in mu=GM! overloaded notation
    for i = 1:N
        mu_vec(i) = 0;
        for j = 1:N
            mu_vec(i) = mu_vec(i) + B_vec(j)*I(j,i);
        end
    end
    p0_vec = zeros(N,1);
    q0_vec = zeros(N,1);
    % eq 7.149, 7.150
    for i = 1:N
        p0_vec(i) = mu_vec(i)/(B_scalar-f(i))*sin(gamma(i));
        q0_vec(i) = mu_vec(i)/(B_scalar-f(i))*cos(gamma(i));
    end
    p0 = -sum(p0_vec);
    q0 = -sum(q0_vec);
    i0 = asin( sqrt(q0^2+p0^2) );
    Omega0 = atan2(p0,q0);
    Omega0 = mod(Omega0,2*pi);
    warning('on')
end
