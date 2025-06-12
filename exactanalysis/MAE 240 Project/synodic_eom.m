function dXdt = synodic_eom(~, X, nu)
    % Extract elements from the input state vector
    x = X(1);
    y = X(2);
    z = X(3);
    vx = X(4);
    vy = X(5);
    vz = X(6);
    
    % Compute squared distances to each primary
    r1_squared = (x + nu)^2 + y^2 + z^2;
    r2_squared = (x - 1 + nu)^2 + y^2 + z^2;
    r1_cu = r1_squared^(3/2);
    r2_cu = r2_squared^(3/2);
    
    % Partial derivatives of pseudo‐potential U (from lecture)
    Ux = x - (1 - nu) * (x + nu) / r1_cu - nu * (x - 1 + nu) / r2_cu;
    Uy = y - (1 - nu) * y / r1_cu - nu * y / r2_cu;
    Uz = -(1 - nu) * z / r1_cu - nu * z / r2_cu;
    
    % Equations of motion in the synodic frame
    d2dx2 = 2 * vy + Ux;
    d2dy2 = -2 * vx + Uy;
    d2dz2 = Uz;
    
    % Return the derivatives of the state vector
    dXdt = [vx; vy; vz; d2dx2; d2dy2; d2dz2];
end