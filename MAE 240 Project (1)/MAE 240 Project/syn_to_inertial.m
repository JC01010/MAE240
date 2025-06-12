function [rv_I] = syn_to_inertial(rv_S , RV_R_barycenter_unitless , const)

% Adopted Professor Rosengren's code for transforming from the synodic
% frame to the inertial frame

%{
MODIFIED FILENAME: syn_to_inertial

DESCRIPTION: This function convert a test particle state vector (position and velocity)
in the synodic rotating frame (dimensionless) to the inertial (ecliptic) ECLIPJ2000 frame

INPUTS:
    rv_S        = the secondary's position and velocity vectors in the inertial frame [ km | km/s ]
    RV_R_barycenter_unitless         
                = position and velocity vectors defined by the rotating frame of the primary and the secondary 
                  with unit length being the distance between the seondary and the primary
    const       = gravitational parameter of primary and secondary bodies [ km^3/s^2 ]

OUTPUTS:
    rv_I        = the target body's position and velocity vector in the inertia frame [ km | km/s ]

REFERENCES:
    https://help.agi.com/stk/index.htm#gator/eq-rlp.htm?Highlight=rotating%20frame
%}

% rotating frame unitless to rotating frame with unit
lstar           = norm( rv_S(1:3) ); % characteristic length
mustar          = const.primary.mu + const.secondary.mu; % G*(characteristic u)
tstar           = sqrt( lstar^3 / mustar );
R_R_barycenter_unitless     = RV_R_barycenter_unitless(1:3); % unitless R in rotating barycenter
V_R_barycenter_unitless     = RV_R_barycenter_unitless(4:6); % unitless V in rotating barycenter
R_R_barycenter  = R_R_barycenter_unitless * lstar; % unit in km R in rotating barycenter
V_R_barycenter  = V_R_barycenter_unitless * lstar / tstar; % unit in km/s V in rotating barycenter
RV_R_barycenter = [ R_R_barycenter ; V_R_barycenter ];

% Take ephemeris of Jupiter position and calculate circular orbit velocity and create DCM
r_S             = rv_S(1:3); %km
v_S             = rv_S(4:6); %km

% r3bp requires circle orbit
omega           = sqrt( const.primary.mu / norm( r_S )^3);
hhat            = cross( r_S , v_S ) / norm( cross( r_S , v_S ) ); % direction of out of plane
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
QdQQ_matrix     = [ Q , zeros(3) ; dQ_matrix , Q ];
QdQQ_matrix_inverse = QdQQ_matrix^-1;

% rotating barycenter frame -> shift barycenter to sun-centered inertial rv -> J2000 sun-centered inertial rv
rv_I = QdQQ_matrix_inverse * RV_R_barycenter + [ R0 ; V0 ];

end

