clear; clc; close all;

% Input mass ratio nu = m_Earth / (m_Sun + m_Earth) from JPL
nu = 3.054200000000000E-6; 

% Initial nondimensionalized position vector from JPL [ km ]
r0_syn = [1.0112931844584114E+0; 2.7024239406843129E-23; 9.5524288608004041E-6];
% Initial nondimensionalized velocity vector from JPL [ km/s ]
v0_syn = [1.1603448296085934E-15; -8.9639462449219390E-3; 3.7508530364328665E-18];
initial_condition = [r0_syn; v0_syn];   % First three -- Position, last three -- Velocity

% Define the nondimensionalized time span for the ode45 integrator
T = 1*2*pi;
tspan  = [0, T];

% Define tolerances for the ode45 integrator
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Integrate with ode45
[t, final] = ode45(@synodic_eom, tspan, initial_condition, options, nu);

% Extract position and velocity components after the integration
x = final(:, 1);
y = final(:, 2);
z = final(:, 3);
vx = final(:, 4);
vy = final(:, 5);
vz = final(:, 6);

% Compute Jacobi constant at each time step
r1 = sqrt( (x + nu).^2 + y.^2 + z.^2 );
r2 = sqrt( (x - 1 + nu).^2 + y.^2 + z.^2 );

% Pseudo‐potential U(x,y,z) in synodic frame (from lecture)
U = 0.5 * ( x.^2 + y.^2 ) + (1 - nu) ./ r1 + nu ./ r2;

% Jacobi constant: J = 2U - (dx^2 + dy^2 + dz^2)
J = 2 * U - ( vx.^2 + vy.^2 + vz.^2 );
fprintf('The Jacobi constant is conserved and results in:\n');
disp(J);

% Plot the XY-plane projection to show the original halo orbit
figure('Position',[200 200 600 600]);
plot(x, y, 'Color',[0.85 0.33 0.10], 'LineWidth',1.2);
hold on;
plot(-nu,  0, 'ro', 'Color',[0.93 0.69 0.13], 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'LineWidth',1.5);
plot(1-nu, 0, 'bo', 'Color',[0 0.45 0.74], 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'LineWidth',1.5);
grid on;
axis equal;
xlabel('x (unitless)');
ylabel('y (unitless)');
title('Original Halo Orbit: XY-plane Projection');
legend('CR3BP trajectory','Sun (primary)','Earth (secondary)', 'Location','best');
hold off;

% Plot the Jacobi constant vs time to confirm that it is a conserved quantity
figure('Position',[850 200 600 400]);
plot(t, J, 'LineWidth',1.3, 'Color',[0.15 0.40 0.68]);
grid on;
xlabel('t (unitless)');
ylabel('Jacobi Constant J');
title('Jacobi Constant vs. Time');
ylim([min(J)*0.9999, max(J)*1.0001]);  % zoom in on variation


%%%%%%%%%%%%%%%%%%%%%%%%%%% SYNODIC TO INERTIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off

% Clear any existing kernels
cspice_kclear()

% Set base path to your desktop
desktop_path = fullfile('/Users', 'shengjingtong', 'Desktop');

% Add MICE paths (assuming standard SPICE installation structure)
addpath(fullfile(desktop_path, 'mice3', 'lib'));
addpath(fullfile(desktop_path, 'mice3', 'src', 'mice'));

% Verify the paths were added
which cspice_furnsh

% Load required kernels
cspice_furnsh(fullfile(desktop_path, 'mice3', 'de430.bsp'));
cspice_furnsh(fullfile(desktop_path, 'mice3', 'pck00010.tpc'));
cspice_furnsh(fullfile(desktop_path, 'mice3', 'gm_de431.tpc'));
cspice_furnsh(fullfile(desktop_path, 'mice3', 'naif0011.tls'));

% Sun & Earth gravitational parameters [ km^3/s^2 ]
const.sun.mu = 1.32712e11;
const.earth.mu = 3.986004354360959e5;

% Define primary and secondary constants
const.primary.mu    = const.sun.mu;
const.secondary.mu  = const.earth.mu;
const.sys.mu        = const.secondary.mu / ( const.secondary.mu + const.primary.mu );

% Dimensionalize length and time
lstar = 149e06;  % km
tstar = sqrt((lstar)^3 / (const.primary.mu + const.secondary.mu));  % seconds
jd_ref = 2451545.0; % Reference Julian date (J2000.0)
const.et0 = cspice_str2et('2000 JAN 01 12:00:00 TDB');  % Reference epoch in ET

% Convert Synodic to Inertial Frame
inertial = zeros(length(t), 6);
state_earth = zeros(length(t), 6);

for k = 1:length(t)
    % Convert non-dimensional time to ET (ephemeris time in seconds past J2000)
    et_current = const.et0 + t(k) * tstar;
    
    % Get Earth state in Sun-centered ecliptic frame
    [state_earth_k, ~] = cspice_spkezr('Earth', et_current, 'ECLIPJ2000', 'NONE', 'Sun');
    state_earth(k,:) = state_earth_k';
    
    % Convert synodic to inertial frame (dimensional)
    inertial(k,:) = syn_to_inertial(state_earth(k,:)', final(k,:)', const)';
end

% Initial conditions for inertial integration at the beginning
r_jwst0 = inertial(1, 1:3)';
v_jwst0 = inertial(1, 4:6)';

% Set up the integration conditions
inertial_condition = [r_jwst0; v_jwst0];

% Set up the time for the integrator
T2 = 10 * 365 * 86400;  % seconds (1 year)
tspan1 = [0, T2];

% Define tolerances for the later ode45 integrator
option = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Integrate with the ode45 integrator
[time, state] = ode45(@inertial_eom, tspan1, inertial_condition, option);

% Prepare Earth state vectors for transformation back to synodic
% Create time vector
et_last = const.et0 + time;  % ET times for integration period

% Get Earth states during integration period
state_earth_new = zeros(length(time), 6);
for k = 1:length(time)
    [state_earth_k, ~] = cspice_spkezr('Earth', et_last(k), 'ECLIPJ2000', 'NONE', 'Sun');
    state_earth_new(k,:) = state_earth_k';
end

% Transform inertial states to synodic barycenter frame
syn = zeros(size(state));
for k = 1:length(time)
    syn(k,:) = inertial_to_syn(state_earth_new(k,:)', state(k,:)', const)';
end

% Plot the XY-plane projection to show the transformed halo orbit in the
% synodic frame
figure('Position',[200 200 600 600]);
plot(x, y, 'Color',[0.85 0.33 0.10], 'LineWidth',1.2);
hold on;
plot(syn(:, 1), syn(:, 2), 'Color',[0.85 0.33 0.10], 'LineWidth',1.2);
hold on;
plot(-nu,  0, 'ro', 'Color',[0.93 0.69 0.13], 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'LineWidth',1.5);
plot(1-nu, 0, 'bo', 'Color',[0 0.45 0.74], 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'LineWidth',1.5);
grid on;
axis equal;
xlabel('x_{new} (unitless)');
ylabel('y_{new} (unitless)');
title('Transformed Halo Orbit: XY-plane Projection');
hold off;

% Compute the new Jacobi constant at each time step
r3 = sqrt( (syn(:, 1) + nu).^2 + syn(:, 2).^2 + syn(:, 3).^2 );
r4 = sqrt( (syn(:, 1) - 1 + nu).^2 + syn(:, 2).^2 + syn(:, 3).^2 );

% Pseudo‐potential U(x,y,z) in synodic frame (from lecture)
Unew = 0.5 * ( syn(:, 1).^2 + syn(:, 2).^2 ) + (1 - nu) ./ r3 + nu ./ r4;

% Jacobi constant: J = 2U - (dx^2 + dy^2 + dz^2)
Jnew = 2 * Unew - ( syn(:, 4).^2 + syn(:, 5).^2 + syn(:, 6).^2 );

% Plot the new Jacobi constant vs time to confirm that it is a conserved quantity
figure('Position',[850 200 600 400]);
plot(time/tstar, Jnew, 'LineWidth',1.3, 'Color',[0.15 0.40 0.68]);
grid on;
xlabel('t (unitless)');
ylabel('New Jacobi Constant J');
title('Updated Jacobi Constant vs. Time');
ylim([min(Jnew)*0.9999, max(Jnew)*1.0001]);  % zoom in on variation

