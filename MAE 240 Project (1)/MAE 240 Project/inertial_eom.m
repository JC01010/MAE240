function dstatedt = inertial_eom(t, state)
    
    % Input Sun and Earth's gravitational paramaters [ km^3/s^2 ]
    mu_s = 1.32712e11;
    mu_e = 3.986004354360959e5;
    
    % Extract elements from the input state vector
    r = state(1:3);
    v = state(4:6);

    % Make sure the "cspice_str2et" is called only once
    persistent t0
    if isempty(t0)
        t0 = cspice_str2et('2000 JAN 1 12:00:00');  % J2000 epoch
    end
    tnew = t0 + t;
    
    % Compute general two body acceleration
    a_2bp = -mu_s / norm(r)^3 * r;
    
    % Get Earth's state (position) from SPICE downloaded from NASA
    [r_earth, ~] = cspice_spkpos('Earth', tnew, 'ECLIPJ2000', 'None', 'Sun');
    r_e = r_earth;

    % Compute the bodycentric frame 3BP of Earth's perturbation
    accel_term = ((r - r_e) / norm(r - r_e)^3) + r_e / norm(r_e)^3;
    a_pert = -mu_e * accel_term;

    % Add in the SRP
    albedo = 0.9; % Albedo of the JWST due to that it is made of aluminum
    % The area of the JWST windshield is as big as a tennis court from JWST
    A_m = 260.87 / 6500;
    P_phi = 4.56e-06; % The solar radiation constant
    beta = (1 + albedo) * A_m * P_phi;
    d = r * 1000; % km to m
    d_hat = d / norm(d);
    a_srp = (beta / norm(d)^2 * d_hat) / 1000; % m back to km

    % Compute the total acceleration
    a_total = a_2bp + a_pert + a_srp;
    
    % Return the derivatives of the state vector
    dstatedt = [v; a_total];
end