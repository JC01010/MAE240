function [RV_R_barycenter_unitless] = inertial_to_syn(rv_S , rv_I , const)

% Adopted Professor Rosengren's code for transforming from the inertial
% frame to the synodic frame

%{
MODIFIED FILENAME: inertial_to_syn

DESCRIPTION: This function convert a test particle state vector (position and velocity)
in inertial (ecliptic) ECLIPJ2000 frame to the synodic rotating frame (dimensionless)

INPUTS:
    rv_S        = the secondary's position and velocity vectors in the inertial frame [ km | km/s ]
    rv_I        = the target body's position and velocity vector in the inertia frame [ km | km/s ]
    const       = gravitational parameter of primary and secondary bodies [ km^3/s^2 ]

OUTPUTS:
    RV_R_barycenter_unitless         
                = position and velocity vectors defined by the rotating frame of the primary and the secondary 
                  with unit length being the distance between the seondary and the primary

REFERENCES:
    https://help.agi.com/stk/index.htm#gator/eq-rlp.htm?Highlight=rotating%20frame
%}

% Ecliptic J2000 primary-centered Inertial frame r and v in column vector
rv_I            = rv_I(:);

% From ephemeris of secondary position, calculate circular orbit velocity and create DCM
rv_S            = rv_S(:);
r_S             = rv_S(1:3); 
v_S             = rv_S(4:6); 

% r3bp requires circle orbit 
omega           = sqrt( const.primary.mu / norm( r_S )^3 );
hhat            = cross(r_S , v_S) / norm( cross( r_S , v_S ) ); % direction of out of plane
omega_vector    = omega * hhat;
v_S_omega       = cross( omega_vector , r_S ); % circular assumption velocity vector

% circular assumption 
v_S_cir         = v_S_omega; % use circular assumption for velocity vector

% instanteneous rotating axes
xhat            = r_S / norm( r_S );
zhat            = cross( r_S , v_S_cir ) / norm( cross( r_S , v_S_cir ) );
yhat            = cross( zhat , xhat );

% position transfrom related
% rotating frame to inertial frame DCM (right to left)
I_DCM_R         = [ xhat , yhat , zhat ];

% inertial frame to rotating frame DCM (Just the transpose of previous one)
R_DCM_I         = I_DCM_R';

% Another way express dxhat dyhat dzhat from AGI page
% benefit no need to mess with angular velocity and its direction
ratio_u         = const.secondary.mu / ( const.primary.mu + const.secondary.mu );
dxhat           = v_S_cir / norm( r_S ) - r_S * ( r_S' * v_S_cir ) / norm( r_S )^3;
dyhat           = cross( cross( r_S , v_S_cir ) , v_S_cir ) / ...
                  ( norm( r_S ) * norm( cross( r_S , v_S_cir ) ) ) - ...
                  ( cross( r_S , v_S_cir ) / ( norm( r_S )^3 * ...
                  norm( cross( r_S , v_S_cir ) ) ) )' * ( cross( cross( r_S , v_S_cir ) , r_S ) );
dzhat           = [ 0 ; 0 ; 0 ];
dQ_matrix       = [ dxhat , dyhat , dzhat ]';
Q               = R_DCM_I;
R0              = ratio_u * r_S;
V0              = ratio_u * v_S_cir;
QdQQ_matrix     = [ Q , zeros(3) ; dQ_matrix,Q ];

% primary-centered inertial rv -> shift primary-centered inertial rv to barycenter -> rotating barycenter frame
r_I             = rv_I(1:3);
v_I             = rv_I(4:6);
RV_R_barycenter = QdQQ_matrix * [ r_I - R0 ; v_I - V0 ];
% RV_R_barycenter = QdQQ_matrix*[r_I;v_I];
R_R_barycenter  = RV_R_barycenter(1:3);
V_R_barycenter  = RV_R_barycenter(4:6);

% rotating frame with unit to rotating frame unitless
lstar           = norm( rv_S(1:3) ); % characteristic length
mustar          = const.primary.mu + const.secondary.mu; % G*(characteristic u)
tstar           = sqrt( lstar^3 / mustar );

R_R_barycenter_unitless     = R_R_barycenter / lstar;
V_R_barycenter_unitless     = V_R_barycenter * tstar/lstar;
RV_R_barycenter_unitless    = [ R_R_barycenter_unitless ; V_R_barycenter_unitless ];

end

