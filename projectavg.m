close all; clearvars; clc;
% Initial conditions (normalized units)
r0 = [1.0112931844584114, 2.7024239406843129E-23, 9.5524288608004041E-6];
v0 = [1.1603448296085934E-15, -8.9639462449219390E-3, 3.7508530364328665E-18];

% r0 and v0 are row or column vectors
r0 = r0(:);
v0 = v0(:);

% Magnitudes
r_mag = norm(r0);
v_mag = norm(v0);

% Specific angular momentum
h_vec = cross(r0, v0);
h_mag = norm(h_vec);

% Semi-major axis (a)
mu = 1; % normalized units
a = 1 / (2/r_mag - v_mag^2/mu);

% Eccentricity vector and magnitude (e0)
e_vec = (cross(v0, h_vec) / mu) - (r0 / r_mag);
e0 = norm(e_vec);

% Argument of periapsis (w0)
% For planar orbits, argument of periapsis is the angle from x-axis to e_vec
w0 = atan2(e_vec(2), e_vec(1));

% Normalized mean motion (n)
n = sqrt(mu / a^3);

% Display results
fprintf('Initial semi-major axis a = %.8f\n', a);
fprintf('Initial eccentricity e0 = %.8f\n', e0);
fprintf('Initial argument of periapsis w0 = %.8f rad\n', w0);
fprintf('Normalized mean motion n = %.8f\n', n);

% Orbital elements
% e0 = 0.01; % initial eccentricity %% UPDATE
% w0 = 0;    % initial argument of periapsis (rad) %% UPDATE
% a = 1;     % normalized semi-major axis %% UPDATE
% n = 1;     % normalized mean motion %% UPDATE

% SRP parameters
P0 = 4.56e-6; % N/m^2
Cr = 1.5;     % reflectivity coefficient
A = 25;     % m^2 (JWST)
m = 6500;     % kg (JWST)

a_SRP_SI = P0 * Cr * A / m; % [m/s^2]

AU = 1.496e11; % m
year = 3.15576e7; % s

aa_SRP = a_SRP_SI * year^2 / AU;

d_hat = [1; 0; 0];
e_hat = [cos(w0); sin(w0); 0];

% Time span
T = 100;
tspan = linspace(0, T, 1000);

odefun = @(t, y) averaged_eom(t, y, n, a, a_SRP_SI, d_hat);

y0 = [e0; w0];
[t, y] = ode45(odefun, tspan, y0);

e = y(:,1);
w = y(:,2);

% Compute secular trends at each time step
dRdw = zeros(size(t));
dRde = zeros(size(t));
for i = 1:length(t)
    e_hat_i = [cos(w(i)); sin(w(i)); 0];
    dotprod = dot(d_hat, e_hat_i);
    dRdw(i) = -1.5 * a_SRP_SI * e(i) * derivative_dot_d_e_wrt_w(w(i), d_hat);
    dRde(i) = -1.5 * a_SRP_SI * dotprod;
end

figure;
plot(t, e, 'b-');
xlabel('Time'); ylabel('Eccentricity e');
title('Evolution of Eccentricity under SRP');

figure;
plot(t, w, 'r-');
xlabel('Time'); ylabel('Argument of Periapsis w (rad)');
title('Evolution of Argument of Periapsis under SRP');

figure;
plot(t, dRdw, 'k-');
xlabel('Time'); ylabel('\partialR/\partialw');
title('Secular Trend: \partialR/\partialw');

figure;
plot(t, dRde, 'm-');
xlabel('Time'); ylabel('\partialR/\partiale');
title('Secular Trend: \partialR/\partiale');

Tvals = [10, 100, 500, 1000];
titles = {'T = 10', 'T = 100', 'T = 500', 'T = 1000'};

% Preallocate storage for results
e_all = cell(1,4);
w_all = cell(1,4);
dRdw_all = cell(1,4);
dRde_all = cell(1,4);
t_all = cell(1,4);

for k = 1:4
    T = Tvals(k);
    tspan = linspace(0, T, 1000);

    y0 = [e0; w0];
    [t, y] = ode45(@(t, y) averaged_eom(t, y, n, a, a_SRP_SI, d_hat), tspan, y0);

    e = y(:,1);
    w = y(:,2);

    % Compute secular trends at each time step
    dRdw = zeros(size(t));
    dRde = zeros(size(t));
    for i = 1:length(t)
        e_hat_i = [cos(w(i)); sin(w(i)); 0];
        dotprod = dot(d_hat, e_hat_i);
        dRdw(i) = -1.5 * a_SRP_SI * e(i) * derivative_dot_d_e_wrt_w(w(i), d_hat);
        dRde(i) = -1.5 * a_SRP_SI * dotprod;
    end

    e_all{k} = e;
    w_all{k} = w;
    dRdw_all{k} = dRdw;
    dRde_all{k} = dRde;
    t_all{k} = t;
end

% Plot grouped by secular trend, not by T
figure;
for k = 1:4
    subplot(2,2,k);
    plot(t_all{k}, e_all{k}, 'b-');
    xlabel('Time'); ylabel('Eccentricity e');
    title(['Eccentricity, ' titles{k}]);
end
sgtitle('Eccentricity for Different T');

figure;
for k = 1:4
    subplot(2,2,k);
    plot(t_all{k}, w_all{k}, 'r-');
    xlabel('Time'); ylabel('Argument of Periapsis w (rad)');
    title(['w, ' titles{k}]);
end
sgtitle('Argument of Periapsis for Different T');

figure;
for k = 1:4
    subplot(2,2,k);
    plot(t_all{k}, dRdw_all{k}, 'k-');
    xlabel('Time'); ylabel('\partialR/\partialw');
    title(['\partialR/\partialw, ' titles{k}]);
end
sgtitle('\partialR/\partialw for Different T');

figure;
for k = 1:4
    subplot(2,2,k);
    plot(t_all{k}, dRde_all{k}, 'm-');
    xlabel('Time'); ylabel('\partialR/\partiale');
    title(['\partialR/\partiale, ' titles{k}]);
end
sgtitle('\partialR/\partiale for Different T');

% Averaged equations of motion
function dydt = averaged_eom(~, y, n, a, aa_SRP, d_hat)
    e = y(1);
    w = y(2);
    e_hat = [cos(w); sin(w); 0];
    e_hat = [cos(w); sin(w); 0];
    dotprod = dot(d_hat, e_hat);
    sqrt1me2 = sqrt(1 - e^2);
    % Partial derivatives

    dRdw = -1.5 * aa_SRP * e * derivative_dot_d_e_wrt_w(w, d_hat);
    dRde = -1.5 * aa_SRP * dotprod;
    % Averaged equations
    % Averaged equationsn * a^2 * e) * dRdw;
    dedt = -sqrt1me2 / (n * a^2 * e) * dRdw;
    dwdt = -sqrt1me2 / (n * a^2 * e) * dRde;
    dydt = [dedt; dwdt];
    dydt = [dedt; dwdt];
end
% Derivative of dot(d_hat, e_hat) with respect to w
% Derivative of dot(d_hat, e_hat) with respect to w
function val = derivative_dot_d_e_wrt_w(w, d_hat)
    % d/dw (d̂∙ê) = d̂∙dê/dw = d̂∙[-sin(w); cos(w); 0]
    val = -d_hat(1)*sin(w) + d_hat(2)*cos(w);end