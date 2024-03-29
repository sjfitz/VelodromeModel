function [Track, Info, Edge] = VelodromeModel(Y, R, n, L_L, Opts)
% VelodromeModel creates a velodrome track black-line model that consists of 
% two straights, two circular arc bends and four transition curves between the 
% bends and the straights. The transition curves are based on different Cesaro
% equations where the curvature along the arc length is defined. This provides
% controlling of the centripetal acceleration during cornering and has chosen 
% levels of geometric continuity. 
% 
% The first option is for a generalised clothoid (also known as a Euler spiral 
%   or Cornu spiral). Here the curvature of the transition curve is a 
%   polynomial function of its arc length and the exponent power, n, is 
%   required to be selected. A typical case is for n = 1 which is the standard 
%   Euler spiral transition curve where the curvature increases linearly along 
%   the arc length. This is a G2 continuous curve.
% 
% The second option is for a half-sine wave curvature profile where the 
%   curvature increases from zero to the bend curvature following a sinusoidal
%   path. This is a G3 continuous curve. This is also called the 'cosine' curve.
% 
% The third option is for a general polynomial curvature profile where the 
%   curvature increases from zero to the bend curvature following a polynomial
%   path. This can be anything from a G2 to G(X) path and the polynomial order
%   depends on the chosen geometric continuity. Higher levels of continuity
%   require longer transition lengths, providing an upper limit to the 
%   maximum achievable continuity. 
% 
% The measurable features that define the track are:
%   Y       The half-span between the two straights.
%   R       The turn radius at the bend apex.
%   L_L     The lap length. This is generally a known, fixed value. 
%   If a clothoid is chosen
%       n   The power of the change in curvature with length. Usually 1 (linear)
% 
% The ratio of inputs Y/R has a limited range of feasible solutions dependent on
% L_L, R, and the curvature function that are checked before calculations begin. 
% 
% Requires Matlab r2019b or later. 
% 
% Syntax
%   [Track, Info, Edge] = VelodromeModel(Y, R)
%   ... = VelodromeModel(Y, R, 1)
%   ... = VelodromeModel(Y, R, 'sine')
%   ... = VelodromeModel(Y, R, 'G3')
%   ... = VelodromeModel(Y, R, 'G4')
%   ... = VelodromeModel(Y, R, 'G5')
%   ... = VelodromeModel(Y, R, n, L_L)
%   ... = VelodromeModel(..., 'Bank',[<value>, <value>])
%   ... = VelodromeModel(..., 'Width',<value>)
%   ... = VelodromeModel(..., 'Resolution',<value>)
%   ... = VelodromeModel(..., 'FileName',<name>)
% 
% Inputs (Ordered)
%   Y           (1 x 1 double)  [m] Half-span between the two straights.
%   R           (1 x 1 double)  [m] Radius of the circular bend arc. R < Y.
%   n       Different curvature profile options:
%               (1 x 1 double)  #       G2 - Power
%                   n > 0. Default 1.
%               (1 x n char)    'sine'  G3 - Sinusoidal
%               (1 x n char)    'g#'    G# - Polynomial
%                   Provides any chosen level of G continuity (up to practical
%                   limits). Higher G levels requires longer transitions. 
%   L_L         (1 x 1 double)  [m] Lap length. Default 250.
% 
% Inputs (Name-value pairs) (optional) 
%   Bank        (1 x 2 double)  [deg] Create a simple sinusoidal bank angle:
%               The minimum and maximum values of the bank angle. 
%                   BankAngle = (Max - Min)/2*sin(2*t - pi/2) + (Max + Min)/2
%   Width       (1 x 1 double)  [m] The track width. 
%   Resolution  (1 x 1 double)  [m] Resolution of output data points. 
%                   Default 1. (Every 1 m along the datum line). 
%   FileName    (1 x n char)    Path/filename/extension to save the output data.
%                   File type based on extension: .txt, .dat, .csv, .xls, .xlsx
% 
% Outputs 
%   Track     	(m x 9+ table) A table of the data vs track distance.
%                       m = L_L/Resolution + 1. 
%               Variables:
%                   Lap         [m]     Lap distance from the pursuit line
%                   X           [m]     x-coordinate of the datum line
%                   Y           [m]     y-coordinate of the datum line
%                   Z           [m]     z-coordinate of the datum line (0)
%                   k           [m^-1]  Curvature
%                   Radius      [m]     Radius of curvature
%                   dk_ds       [m^-2]  1st derivative of curvature w.r.t. Lap
%                   d2k_ds2     [m^-3]  2nd derivative of curvature w.r.t. Lap
%                   Tangent     [rad]   Tangential angle
%               If 'Bank' is included
%                   BankAngle   [deg]   Banking angle
%               If 'Bank' & 'Width' are included
%                   X_Top       [m]     x-coordinate of the track top
%                   Y_Top       [m]     y-coordinate of the track top
%                   Z_Top       [m]     z-coordinate of the track top
%               Track also has two structures 'Info' and 'Edge' in the table
%               custom properties that record calculation details. 
%   Info        (1 x 1 struct) Various track information.
%   Edge        (1 x 1 struct) The (x, y) location of the track change points. 
% 
% Example usage:
%   Track = VelodromeModel(23, 22);
%   Track = VelodromeModel(23, 22, 'sine', 250, 'Bank',[13, 43], 'Width',7.5);
%   Track = VelodromeModel(23, 22, 'G4');
%   Track = VelodromeModel(23, 22,    1,   250, 'Resolution',0.25);
%   Track = VelodromeModel(23, 22,    1,   250, 'FileName','TrackData.csv');
%   figure; stackedplot(Track, 'XVariable','Lap');
%   figure; plot(Track.X, Track.Y); axis equal;
%   figure; plot(Track.Lap, Track.k); 
% 
% Corresponding article:
%   'Impact of transition design on the accuracy of velodrome models' 
%   Sports Engineering, 24(23), https://doi.org/10.1007/s12283-021-00360-3.
%   Fitzgerald, S.(1), Kelso, R.(1), Grimshaw, P.(1,2), Warr, A.(3) 
%   1 School of Mechanical Engineering, The University of Adelaide, Australia
%   2 College of Health and Life Sciences, Hamad Bin Khalifa University, Qatar
%   3 Adelaide, Australia
%   Contact: shaun.fitzgerald@adelaide.edu.au 
% Derivation of 'G# - Polynomial' option in PhD Thesis: 
%   'Understanding the Aerodynamic Environment in Track Cycling'
%   Shaun Fitzgerald, The University of Adelaide. 
% 
% Maintained at: https://github.com/sjfitz/VelodromeModel

%% Inputs 
arguments
    Y               (1,1) double {mustBePositive}
    R               (1,1) double {mustBePositive}
    n               (1,:) {double, char} = 1
    L_L             (1,1) double {mustBePositive} = 250
    Opts.Bank       (1,2) double {mustBeNonnegative} = [0, 0]
    Opts.Width      (1,1) double {mustBeNonnegative} = 0
    Opts.Resolution (1,1) double {mustBePositive} = 1
    Opts.FileName   (1,:) {char, string}
end

nDataP = 500; % [#] Number of data points for internal calculations

% Basic bounds
assert(Y < L_L/(2*pi), 'The half-span Y must be < L_L/(2*pi).')
assert(R < Y, 'The radius R must be < Y.')
assert(Opts.Resolution <= L_L/25, 'The resolution must be << L_Lap.')

%% Transition curve calculations 
Input_n = n;
if isnumeric(n) 
    %% Curvature profile: Power equation (17) 
    assert(isscalar(n) && n > 0, 'The exponent n must be scalar and > 0.')
    
    Style = sprintf('Power, n: %g', n);
    Continuity = 'G2';
    
    % Functions
    K  = @(v, L_T, R) v.^n/(R*L_T^n);                           % Curvature   
    Kd = @(v, L_T, R) n*v.^(n-1)/(R*L_T^n);                     % 1st derivative 
    K2d= @(v, L_T, R) n*(n-1)*v.^(n-2)/(R*L_T^n);               % 2nd derivative 
    Ki = @(u, L_T, R) u.^(n+1)/(R*L_T^n*(n+1));                 % Integral
    IC = @(t, L_T, R) integral(@(u) cos(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    IS = @(t, L_T, R) integral(@(u) sin(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt((n+2)*(n+1)*R*(Y - R));
    
    % Saving the curvature function as a string
    CurvFunc = sprintf('K(t) = 1/(R*L_T^%g)*t^%g', n, n);
    if n == 1
        CurvFunc_Latex = 'K(\tau) = \frac{\tau}{RL_T}';
    else
        CurvFunc_Latex = [...
            'K(\tau) = \frac{1}{R}\left(\frac{\tau}{L_T}\right)^{', ...
            num2str(n, '%g'), '}'];
    end
    
elseif strcmpi(n(1), 's')
    %% Curvature profile: Sinusoidal (18) 
    Style = 'Sinusoidal';
    Continuity = 'G3';
    
    % Functions
    K  = @(v, L_T, R) (sin(pi/L_T*v - pi/2) + 1)/2/R;           % Curvature 
    Kd = @(v, L_T, R)  cos(pi/L_T*v - pi/2)*pi/(2*R*L_T);       % 1st derivative 
    K2d= @(v, L_T, R) -sin(pi/L_T*v - pi/2)*pi^2/(2*R*L_T^2);   % 2nd derivative 
    Ki = @(u, L_T, R) (-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R;   % Integral
    IC = @(t, L_T, R) integral(@(u) cos((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    IS = @(t, L_T, R) integral(@(u) sin((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    n  = 1;
    
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt(4*pi/(pi - 2)*R*(Y - R));
    
    % Saving the curvature function as a string
    CurvFunc = 'K(t) = 1/(2*R)*(sin(pi/L_T*t - pi/2) + 1)';
    CurvFunc_Latex = ...
        'K(\tau) = \frac{1}{2R} \left(\sin\left(\frac{\pi}{L_T}\tau - \frac{\pi}{2}\right) + 1 \right)';
    
elseif strcmpi(n(1), 'p') || strcmpi(n(1), 'g')
    %% Curvature profile: pth-degree polynomial 
    g = str2double(n(2:end));       % Geometric continuity 
    if g < 2
        error('General polynomials only works for continuity 2 or greater.')
    elseif isnan(g)
        error('Incorrect input for n.')
    end
    p = 2*g - 3;                    % Polynomial order
    Style = sprintf('%gth degree polynomial', p);
    Style = strrep(Style, '1th', '1st');
    Style = strrep(Style, '3th', '3rd');
    Style = strrep(Style, '11st', '11th');
    Style = strrep(Style, '13rd', '13th');
    Continuity = sprintf('G%g', g);
    
    Arr_N = (0:g - 2);
    Arr_A = factorial(p - Arr_N)./factorial(p - Arr_N - Arr_N');
    Arr_B = [1; zeros(g - 2, 1)];
    Ci = Arr_A\Arr_B;
    Ci_All = Arr_A.*Ci';
    
    % Special conditions for the early functions where the coefficient array
    % is smaller than the number of derivatives to be output. 
    if g == 2
        Ci_All(2,:) = Ci_All(1,:);
        Ci_All(3,:) = 0;
    elseif g == 3
        Ci_All(3,:) = Ci_All(2,:).*[2, 1];
    end
    
    % Functions
    K  = @(v, L_T, R) 1/(R*L_T^(g-1))*sum(Ci_All(1,:)./L_T.^(g-Arr_N-2).*repmat(v,1,g-1).^(p-Arr_N),   2, 'omitnan');
    Kd = @(v, L_T, R) 1/(R*L_T^(g-1))*sum(Ci_All(2,:)./L_T.^(g-Arr_N-2).*repmat(v,1,g-1).^(p-Arr_N-1), 2, 'omitnan');
    K2d= @(v, L_T, R) 1/(R*L_T^(g-1))*sum(Ci_All(3,:)./L_T.^(g-Arr_N-2).*repmat(v,1,g-1).^(p-Arr_N-2), 2, 'omitnan');
    Ki = @(u, L_T, R) 1/(R*L_T^(g-1))*sum(Ci_All(1,:)./L_T.^(g-Arr_N-2).*repmat(u,1,g-1).^(p-Arr_N+1)./(p-Arr_N+1), 2);
    IC = @(t, L_T, R) integral(@(u) cos(1/(R*L_T^(g-1))*sum(Ci_All(1,:)./L_T.^(g-Arr_N-2).*repmat(u,1,g-1).^(p-Arr_N+1)./(p-Arr_N+1), 2, 'omitnan')), 0, t, 'ArrayValued',true);
    IS = @(t, L_T, R) integral(@(u) sin(1/(R*L_T^(g-1))*sum(Ci_All(1,:)./L_T.^(g-Arr_N-2).*repmat(u,1,g-1).^(p-Arr_N+1)./(p-Arr_N+1), 2, 'omitnan')), 0, t, 'ArrayValued',true);
    n  = 1;
    
    % The initial value is found with the first-degree Maclaurin polynomial
    alpha = sum(Ci'./((p - Arr_N + 1).*(p - Arr_N + 2)));
    A0 = sqrt(1/alpha*R*(Y - R));
    
    % Saving the curvature function as a string
    CurvFunc = sprintf('K(t) = 1/(R*L_T^%g)*(', g-1);
    for ii = 1:g-1
        CurvFunc = sprintf('%s%.0f/L_T^%g*t^%g + ', ...
            CurvFunc, Ci_All(1,ii), g-ii-1, p-ii+1);
    end
    CurvFunc = sprintf('%s)', CurvFunc(1:end-3));
    
    % Saving the curvature function as a string formatted for latex
    CurvFunc_Latex = sprintf('%s%g%s', ...
        'K(\tau) = \frac{1}{RL_T^', g-1, '} \left( ');
    for ii = 1:g-2
        CurvFunc_Latex = sprintf('%s%s%g%s%g%s{%g} + ', ...
            CurvFunc_Latex, '\frac{', Ci_All(1, ii), '}{L_T^', ...
            g-ii-1, '}\tau^', 2*g-ii-2);
    end
    ii = ii + 1;
    CurvFunc_Latex = sprintf('%s%g%s{%g}', ...
        CurvFunc_Latex, Ci_All(1, ii), '\tau^', 2*g-ii-2);
    CurvFunc_Latex = sprintf('%s%s', CurvFunc_Latex, ' \right)');
    CurvFunc_Latex = strrep(CurvFunc_Latex, '^{}', '');
    CurvFunc_Latex = strrep(CurvFunc_Latex, 'L_T^1}', 'L_T}');
    CurvFunc_Latex = strrep(CurvFunc_Latex, '+ \frac{-', '- \frac{');
    CurvFunc_Latex = strrep(CurvFunc_Latex, '\left( \frac{-', '\left( -\frac{');
    
else
    error('Unknown model type: ''%s''.', n)
    
end

%% Calculations 
% Primary bounds for Y/R
if R <= L_L/(2*pi*(n+1))
    A_Max = pi*R/2*(n+1);                   % [-]    (15) (psi_1 < pi/2)
else
    A_Max = (L_L - 2*pi*R)*(n+1)/(4*n);     % [-]    (16) (L_S > 0)
end
YonR_Max = K(A_Max, A_Max, R)*IS(A_Max, A_Max, R) + cos(Ki(A_Max, A_Max, R));
assert(1 < Y/R && Y/R < YonR_Max, sprintf(...
    'Y/R (%.4f) must be in: (1 < Y/R < %.4f) for ''%s'' curvature with R = %.4g.',...
    Y/R, YonR_Max, Style, R))

% ( 9) Solving for A (tau_1, L_T) with the Newton�Raphson method.
Er = 1;
cc = 0;
while abs(Er) > 1e-10
    A1 = A0 - (IS(A0, A0, R) + R*cos(A0/R/(n+1)) - Y)/(sin(A0/R/(n+1))/(n+1));
    Er = A1 - A0;
    A0 = A1;
    cc = cc + 1;
    if cc > 2000, error('Calculation of A not converged'); end
end
L_T = A1;                                   % [m]    (10) Transition length

% Single-value results  
X       = IC(L_T,L_T,R) - R*sin(Ki(L_T,L_T,R)); % [m]( 9) Bend centre X-coord
psi_1   = Ki(L_T,L_T,R);                    % [rad]  ( 3) Bend start tangent 
theta   = pi/2 - psi_1;                     % [rad]       Bend open angle
L_B     = theta*R;                          % [m]    (11) Bend length
L_S     = L_L/4 - L_T - L_B;                % [m]    (12) Straight length

% Using uneven coordinates so that they are biased towards the joins
rho = linspace(0, pi, nDataP)';
ArrayA = (1 - cos(rho))/2; 
ArrayB = (1 - cos(rho)) - 1;
% ArrayA = linspace( 0, 1, nDataP)';          % Linear method
% ArrayB = linspace(-1, 1, nDataP)';

% Per parametric distance t
t = L_T*ArrayA;                             % [m]    Eqn parameter
for ii = 1:nDataP
    x(ii,1) = IC(t(ii), L_T, R);            % [m]    ( 1) x coordinate
    y(ii,1) = IS(t(ii), L_T, R);            % [m]    ( 1) y coordinate
end
s       = t;                                % [m]    ( 2) Arc length
psi     = Ki( t, L_T, R);                   % [rad]  ( 3) Tangential angle
kappa   = K(  t, L_T, R);                   % [m^-1] ( 4) Curvature
dk_ds   = Kd( t, L_T, R);                   % [m^-2] ( 5) Curvature 1st derivative
d2k_ds2 = K2d(t, L_T, R);                   % [m^-3]      Curvature 2nd derivative
d2k_ds2(isnan(d2k_ds2)) = 0;

% The lap position of each change between track segments 
Transition = cumsum([0, L_S, L_T, 2*L_B, L_T, 2*L_S, L_T, 2*L_B, L_T, L_S]);

%% The (x, y) coordinates for each track segment
% Moving the origin to the centre of the velodrome
xT  = x + L_S;
yT  = y - Y;
xBc = X + L_S;
yBc = 0;

% Edge points of the straights 
Edge.Str_Edge(1,:) = [ L_S, -Y];
Edge.Str_Edge(2,:) = [ L_S,  Y]; 
Edge.Str_Edge(3,:) = [-L_S,  Y]; 
Edge.Str_Edge(4,:) = [-L_S, -Y]; 

% Edge points of the circular bend arcs
Edge.Bnd_Edge(1,:) = [ xBc + R*cos(theta), yBc - R*sin(theta)];
Edge.Bnd_Edge(2,:) = [ xBc + R*cos(theta), yBc + R*sin(theta)];
Edge.Bnd_Edge(3,:) = [-xBc - R*cos(theta), yBc + R*sin(theta)];
Edge.Bnd_Edge(4,:) = [-xBc - R*cos(theta), yBc - R*sin(theta)];

% Centre of the circular bend arc
Edge.Bnd_Centre(1,:) = [ xBc, yBc];
Edge.Bnd_Centre(2,:) = [-xBc, yBc];

% Transition curves
XY.Trn1 = [ xT,  yT];
XY.Trn2 = flipud([ xT, -yT]);
XY.Trn3 = [-xT, -yT];
XY.Trn4 = flipud([-xT,  yT]);

% Straights 
s_S = L_S*ArrayA;                               % [m] Arc length
XY.Str1 = [ L_S*ArrayA,          -Y*ones(nDataP, 1)];
XY.Str2 = [-L_S*ArrayB,           Y*ones(nDataP, 1)];
XY.Str3 = [-L_S*flipud(ArrayA),  -Y*ones(nDataP, 1)];

% Circular bends 
th = 2*theta*ArrayA;                            % [rad] Arc angle
XY.Bnd1 = [ xBc + R*cos(th - theta),        yBc + R*sin(th - theta)];
XY.Bnd2 = [-xBc + R*cos(th - theta + pi),   yBc + R*sin(th - theta + pi)];

%% Lap distance and subsequent parameters
%%%%% Lap distance 
Lap.Str1 = s_S;
Lap.Str2 = 1*L_S + 2*L_T + 2*L_B + 2*s_S;
Lap.Str3 = 3*L_S + 4*L_T + 4*L_B + 1*s_S;

Lap.Bnd1 = 1*L_S + 1*L_T + R*th;
Lap.Bnd2 = 3*L_S + 3*L_T + 2*L_B + R*th;

Lap.Trn1 = L_S + s;
Lap.Trn2 = 1*L_S + 1*L_T + 2*L_B + s;
Lap.Trn3 = 3*L_S + 2*L_T + 2*L_B + s;
Lap.Trn4 = 3*L_S + 3*L_T + 4*L_B + s;

%%%%% Curvature 
Curv.Str1 = zeros(nDataP, 1);
Curv.Str2 = zeros(nDataP, 1);
Curv.Str3 = zeros(nDataP, 1);

Curv.Bnd1 = ones(nDataP, 1)/R;
Curv.Bnd2 = ones(nDataP, 1)/R;

Curv.Trn1 = kappa;
Curv.Trn2 = flipud(kappa);
Curv.Trn3 = kappa;
Curv.Trn4 = flipud(kappa);

%%%%% Curvature 1st derivative
CurvOnDs.Str1 = zeros(nDataP,1);
CurvOnDs.Str2 = zeros(nDataP,1);
CurvOnDs.Str3 = zeros(nDataP,1);

CurvOnDs.Bnd1 = zeros(nDataP,1);
CurvOnDs.Bnd2 = zeros(nDataP,1);

CurvOnDs.Trn1 =  dk_ds;
CurvOnDs.Trn2 = -flipud(dk_ds);
CurvOnDs.Trn3 =  dk_ds;
CurvOnDs.Trn4 = -flipud(dk_ds);

%%%%% Curvature 2nd derivative
CurvOnDs2.Trn1 =  d2k_ds2;
CurvOnDs2.Trn2 =  flipud(d2k_ds2);
CurvOnDs2.Trn3 =  d2k_ds2;
CurvOnDs2.Trn4 =  flipud(d2k_ds2);

%%%%% Tangential angle 
Tangent.Str1 = zeros(nDataP,1);
Tangent.Str2 = zeros(nDataP,1) +   pi;
Tangent.Str3 = zeros(nDataP,1) + 2*pi;

Tangent.Bnd1 = th - theta + 1/2*pi;
Tangent.Bnd2 = th - theta + 3/2*pi;

Tangent.Trn1 = psi;
Tangent.Trn2 = -flipud(psi) + 2*theta + 2*psi_1;
Tangent.Trn3 = psi + pi;
Tangent.Trn4 = -flipud(psi) + 2*theta + 2*psi_1 + pi;

%% Combining the data into one (unequally spaced) lap distance 
Comb.X = [...
    XY.Str1(1:end-1,1); ...
    XY.Trn1(1:end-1,1); ...
    XY.Bnd1(1:end-1,1); ...
    XY.Trn2(1:end-1,1); ...
    XY.Str2(1:end-1,1); ...
    XY.Trn3(1:end-1,1); ...
    XY.Bnd2(1:end-1,1); ...
    XY.Trn4(1:end-1,1); ...
    XY.Str3(1:end-1,1)];

Comb.Y = [...
    XY.Str1(1:end-1,2); ...
    XY.Trn1(1:end-1,2); ...
    XY.Bnd1(1:end-1,2); ...
    XY.Trn2(1:end-1,2); ...
    XY.Str2(1:end-1,2); ...
    XY.Trn3(1:end-1,2); ...
    XY.Bnd2(1:end-1,2); ...
    XY.Trn4(1:end-1,2); ...
    XY.Str3(1:end-1,2)];

Comb.Lap = [...
    Lap.Str1(1:end-1); ...
    Lap.Trn1(1:end-1); ...
    Lap.Bnd1(1:end-1); ...
    Lap.Trn2(1:end-1); ...
    Lap.Str2(1:end-1); ...
    Lap.Trn3(1:end-1); ...
    Lap.Bnd2(1:end-1); ...
    Lap.Trn4(1:end-1); ...
    Lap.Str3(1:end-1)];

Comb.k = [...
    Curv.Str1(1:end-1); ...
    Curv.Trn1(1:end-1); ...
    Curv.Bnd1(1:end-1); ...
    Curv.Trn2(1:end-1); ...
    Curv.Str2(1:end-1); ...
    Curv.Trn3(1:end-1); ...
    Curv.Bnd2(1:end-1); ...
    Curv.Trn4(1:end-1); ...
    Curv.Str3(1:end-1)];

Comb.dk_ds = [...
    CurvOnDs.Str1(1:end-1); ...
    CurvOnDs.Trn1(1:end-1); ...
    CurvOnDs.Bnd1(1:end-1); ...
    CurvOnDs.Trn2(1:end-1); ...
    CurvOnDs.Str2(1:end-1); ...
    CurvOnDs.Trn3(1:end-1); ...
    CurvOnDs.Bnd2(1:end-1); ...
    CurvOnDs.Trn4(1:end-1); ...
    CurvOnDs.Str3(1:end-1)];

Comb.d2k_ds2 = [...
    CurvOnDs.Str1(1:end-1); ...
    CurvOnDs2.Trn1(1:end-1); ...
    CurvOnDs.Bnd1(1:end-1); ...
    CurvOnDs2.Trn2(1:end-1); ...
    CurvOnDs.Str2(1:end-1); ...
    CurvOnDs2.Trn3(1:end-1); ...
    CurvOnDs.Bnd2(1:end-1); ...
    CurvOnDs2.Trn4(1:end-1); ...
    CurvOnDs.Str3(1:end-1)];

Comb.Tangent = [...
    Tangent.Str1(1:end-1); ...
    Tangent.Trn1(1:end-1); ...
    Tangent.Bnd1(1:end-1); ...
    Tangent.Trn2(1:end-1); ...
    Tangent.Str2(1:end-1); ...
    Tangent.Trn3(1:end-1); ...
    Tangent.Bnd2(1:end-1); ...
    Tangent.Trn4(1:end-1); ...
    Tangent.Str3(1:end-1)];

%% Saving the primary information
Info.Style      = Style;
Info.Input_n    = Input_n;
Info.n          = n;
Info.Y          = Y;
Info.X          = X;
Info.R_Bnd      = R;
Info.L_Lap      = L_L;
Info.L_Str      = L_S;
Info.L_Trn      = L_T;
Info.L_Bnd      = L_B;
Info.theta      = theta;
Info.psi_1      = psi_1;
Info.Max_x      = X + L_S + R;
Info.BendCentre = [xBc, yBc];
Info.Transition = Transition;
Info.Continuity = Continuity;
Info.Bank       = Opts.Bank;
Info.Width      = Opts.Width;
Info.Resolution = Opts.Resolution;
Info.CurvFunc   = CurvFunc;
Info.CurvFunc_Latex = CurvFunc_Latex;
if strcmpi(Input_n(1), 'p') || strcmpi(Input_n(1), 'g')
    Info.Ci_All = Ci_All;
else
    Info.Ci_All = nan;
end

%% Creating a consistently spaced table with the data
Track = table;
Track.Lap       = (0:Opts.Resolution:L_L)';
Track.X         = interp1(Comb.Lap, Comb.X,         Track.Lap, 'makima');
Track.Y         = interp1(Comb.Lap, Comb.Y,         Track.Lap, 'makima');
Track.Z         = zeros(height(Track),1);
Track.k         = interp1(Comb.Lap, Comb.k, Track.Lap, 'makima');
Track.Radius    = 1./Track.k;
Track.dk_ds     = interp1(Comb.Lap, Comb.dk_ds,     Track.Lap, 'makima');
Track.d2k_ds2   = interp1(Comb.Lap, Comb.d2k_ds2,   Track.Lap, 'makima');
Track.Tangent   = interp1(Comb.Lap, Comb.Tangent,   Track.Lap, 'makima');

% Adjusting the tangential angle range to (-pi, +pi)
Track.Tangent(round(end/2):end) = Track.Tangent(round(end/2):end) - 2*pi;
% Floating-point error fix
Track.Tangent(abs(Track.Tangent)        < 1e-14) = 0;  
Track.Tangent(abs(Track.Tangent - pi)   < 1e-14) =  pi;
Track.Tangent(abs(Track.Tangent + pi)   < 1e-14) = -pi;

% Setting the track meta data
Track.Properties.VariableUnits = {'m','m','m','m','m^-1','m','m^-2','m^-3','rad'};
Track = addprop(Track, {'Info', 'Edge'}, repmat({'table'},2,1)); 
Track.Properties.CustomProperties.Info = Info;
Track.Properties.CustomProperties.Edge = Edge;

%% Bank Angle 
if ~isequal(Opts.Bank, [0, 0])
    Track.BankAngle = abs(diff(Opts.Bank))/2*sin(4*pi/L_L*Track.Lap-pi/2) + ...
        mean(Opts.Bank);                            % [deg] (19) Bank angle
    Track.Properties.VariableUnits(end) = {'deg'};
    
    % (x, y, z) coordinates of the top of the track
    if Opts.Width ~= 0
        Track.X_Top = Track.X + Opts.Width*cosd(Track.BankAngle).*sin(Track.Tangent);
        Track.Y_Top = Track.Y - Opts.Width*cosd(Track.BankAngle).*cos(Track.Tangent);
        Track.Z_Top = Track.Z + Opts.Width*sind(Track.BankAngle);
        Track.Properties.VariableUnits(end-2:end) = {'m','m','m'};
    end
end

%% Saving the data to a file
if ismember('FileName', fieldnames(Opts))
    writetable(Track, Opts.FileName);
end
