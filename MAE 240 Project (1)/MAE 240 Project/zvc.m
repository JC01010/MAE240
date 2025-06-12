% Adopted Professor Rosengren's code for the Zero Velocity Cruves for CR3BP

%  Earth & Moon gravitational parameters [ km^3/s^2 ]
const.earth.mu = 3.986004354360959e5;          
const.sun.mu = 1.32712e11;         

%  Earth-Sun mass ratio
const.sys.mu    = const.earth.mu / ( const.earth.mu + const.sun.mu );

%  Lagrange equilibrium points 
const.sys.lagrange = lagrange( const.sys.mu );


% ZERO-VELOCITY CURVES BANNER
close all

% Define corresponding pastel background color for each case
pastel_colors = {
    rgb('LightCyan') , ...      % Case 1: Leftmost region
    rgb('MistyRose') , ...      % Case 2: Between 3.1883–3.1722
    rgb('HoneyDew') ,  ...      % Case 3: Between 3.1722–3.0122
    rgb('Lavender') ,  ...      % Case 4: Between 3.0122–2.988
    rgb('PapayaWhip')           % Case 5: Above 2.988
};

% --- Jacobi Integrals ---
[ L1.jacobi , L1.hamil ] = jacobi_C(const.sys.mu, const.sys.lagrange.L1, zeros(3,1));
[ L2.jacobi , L2.hamil ] = jacobi_C(const.sys.mu, const.sys.lagrange.L2, zeros(3,1));
[ L3.jacobi , L3.hamil ] = jacobi_C(const.sys.mu, const.sys.lagrange.L3, zeros(3,1));
[ L4.jacobi , L4.hamil ] = jacobi_C(const.sys.mu, const.sys.lagrange.L4, zeros(3,1));
[ L5.jacobi , L5.hamil ] = jacobi_C(const.sys.mu, const.sys.lagrange.L5, zeros(3,1));

% --- Mesh and ZVC Potential ---
[ x , y ] = meshgrid(-1.5 : 0.005 : 1.5, -1.5 : 0.005 : 1.5);
zvc0 = zeros( size( x ) );
for i = 1:size(x,1)

    for j = 1:size(x,2)

        r = [ x(i,j) ; y(i,j) ; 0 ];
        zvc0(i,j) = -pseudo_U( const.sys.mu , r );

    end
end

% --- Cases ---
hamils = [
    L1.hamil - 0.005;
    L1.hamil + 0.005;
    L2.hamil + 0.005;
    L3.hamil + 0.005;
    L4.hamil + 0.005
];

case_labels = {
    'I : $C > C_1$'
    'II : $C_1 > C > C_2$'
    'III : $C_2 > C > C_3$'
    'IV : $C_3 > C > C_4 = C_5$'
    'V : $C < C_4 = C_5$'
};

% --- Banner Layout ---
figure(10); clf;
set( gcf , 'Color' , 'w' , 'Units' , 'inches' , 'Position' , [ 2 , 6.5 , 9.5 , 2.0 ] );

ha = tight_subplot( 1 , 5 , [ 0.025 0.015 ] , [ 0.05 0.02 ] , [ 0.015 0.015 ] );

for k = 1 : 5
    axes(ha(k)); 
    hold on;

    % Fill subplot background with corresponding pastel
    fill( [ -1.5 1.5 1.5 -1.5 ] , [ -1.5 -1.5 1.5 1.5 ] , pastel_colors{k} , ...
          'EdgeColor' , 'none' );

    % ZVC fill region
    clevel = [ hamils(k) hamils(k) ];
    contourf( x , y , zvc0 , clevel , 'LineWidth' , 1.0 );
    colormap( rgb( 'Gainsboro' ) );

    if k == 5

        % --- Add circular Moon orbit around barycenter (dotted) ---
        theta = linspace( 0 , 2 * pi , 500 );
        r_earth = norm( [ 1 - const.sys.mu , 0 ] );  % Distance from barycenter
        x_circ = r_earth * cos(theta);
        y_circ = r_earth * sin(theta);
        plot( x_circ , y_circ , 'k:' , 'LineWidth' , 1.0 ); 

    end

    % Earth and Moon
    plot( -const.sys.mu , 0 , 'o' , ...
          'MarkerSize' , 8 , 'MarkerFaceColor' , rgb( 'MediumBlue' ) , ...
          'MarkerEdgeColor' , rgb( 'MediumBlue' ) );
    plot( 1 - const.sys.mu, 0 , 'o' , ...
          'MarkerSize' , 2 , 'MarkerFaceColor' , rgb( 'DarkSlateGray' ) , ...
          'MarkerEdgeColor' , rgb( 'DarkSlateGray' ) );

    % Lagrange points in final frame
    if k == 5

        L_pts = { const.sys.lagrange.L1 , const.sys.lagrange.L2 , const.sys.lagrange.L3 , const.sys.lagrange.L4 , const.sys.lagrange.L5 };

        for iL = 1 : numel( L_pts )

            % --- Add location of Lagrange points ---
            plot( L_pts{iL}(1) , L_pts{iL}(2) , 'p' , 'MarkerSize' , 1.0 , ...
                  'MarkerFaceColor' , rgb( 'DarkRed' ) , 'MarkerEdgeColor' , rgb( 'DarkRed' ) );
        end
    end

    axis equal;
    xlim( [ -1.5 1.5 ] ); ylim( [ -1.5 1.5 ] );
    set( gca , 'XColor' , 'none' , 'YColor' , 'none' );
    box on;

    % Frame
    rectangle( 'Position' , [ -1.5 , -1.5 , 3 , 3 ] , 'EdgeColor' , 'k' , 'LineWidth' , 1.25 );

    % Title
    title( case_labels{k} , 'Interpreter' , 'latex' , 'FontSize' , 11 , 'FontWeight' , 'normal' );
end

% --- Export High-Quality PDF ---
    set( gcf , 'Color' , 'w' );
    set( gcf , 'PaperPositionMode' , 'auto' );     
    set( gcf , 'units' , 'inches' , 'NumberTitle' , 'on' , 'Name' , 'Catalog' );
    set( gcf , 'position' , [6,6.75,9.5,5.0] ); 
    set( gcf , 'renderer' , 'Painters' ) %#ok<*FGREN>

% --- Functions ---    
function lagrange = lagrange( mu )

x1_init         = 1 - ( mu / 3 )^( 1/3 );
L1x             = fzero(@(x) equilibrium_function( x , mu ) , x1_init );
x2_init         = 1 + ( mu / 3 )^( 1/3 );
L2x             = fzero(@(x) equilibrium_function( x , mu ) , x2_init );
x3_init         = -1 -( mu / 3 * ( sqrt(2) - 1 ) );
L3x             = fzero(@(x) equilibrium_function( x , mu ) , x3_init );
L45x            = 0.5 - mu;
L4y             = sqrt(3) / 2;
L5y             = -sqrt(3) / 2;

lagrange.L1     = [ L1x , 0 , 0 ]';
lagrange.L2     = [ L2x , 0 , 0 ]';
lagrange.L3     = [ L3x , 0 , 0 ]';
lagrange.L4     = [ L45x , L4y , 0 ]';
lagrange.L5     = [ L45x , L5y , 0 ]'; 

end

function f = equilibrium_function( x , mu )

f = x - ( 1 - mu ) * ( x + mu ) / ( ( abs( x + mu ) )^3 ) - ...
    mu * ( x - 1 + mu ) / ( abs( x - 1 + mu ) )^3;
    
end

function [ jacobi , hamil ] = jacobi_C( mu , r , v )

    r = r(:); v = v(:);
    Ubar = pseudo_U( mu , r );
    jacobi = 2 * Ubar - norm(v)^2;
    hamil = -0.5 * jacobi;

end

function Ubar = pseudo_U( mu , r )

    r1 = [ -mu ; 0 ; 0 ];
    r2 = [ 1.0 - mu ; 0 ; 0 ];
    Ubar = 0.5 * ( r(1)^2 + r(2)^2 ) + ...
           ( 1 - mu ) / norm( r - r1 ) + mu / norm( r - r2 );

end
