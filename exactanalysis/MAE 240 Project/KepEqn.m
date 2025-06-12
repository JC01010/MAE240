function [ f ] = KepEqn( M, e, angletype )
%{
FILENAME: KepEqn

DESCRIPTION: Solves Kepler's equation 

INPUTS:
    M           = mean anomaly at epoch 
    e           = eccentricity
    angletype   = string indicating desired angletype for inputs
                    'deg' = degrees, 'rad' = radians

OUTPUTS:
    f           = true anomaly at epoch

FUNCTIONS/SCRIPTS CALLED:
    none

CALLED BY:
    mean2osc_schaub.m

MODIFICATIONS:
    10-Oct-2010    |    Aaron Rosengren - Original
%}

%  Conversions
d2r = pi/180;

if(strcmp(angletype, 'deg'))
    M = M*d2r;
end

M = mod(M, 2*pi);

%  Compute eccentric anomaly
%  A priori estimate
j = 1;
if -pi < M < 0 || M > pi
    E(j) = M - e;
else
    E(j) = M + e;
end

%  Newton iteration to find eccentric anomaly
%  Algorithm [goal: find E so f = 0]
f_E(j) = E(j) - e*sin(E(j)) - M;
while abs(f_E(j)) > 1e-13
    E(j + 1) = E(j) - f_E(j)/(1 - e*cos(E(j)));
    j = j + 1;
    f_E(j) = E(j) - e*sin(E(j)) - M;
end

%  Converged eccentric anomaly [ rad ]
E = E(j);

%  True anomaly [ rad ]
f = mod(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2)), 2*pi);
% f = mod(mod(f, 2*pi), 2*pi);

if(strcmp(angletype, 'deg'))
    f = f/d2r;
end

end % ---- End Function