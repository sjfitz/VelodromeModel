function Track = VelodromeModel(Y, R, n, L_L, Resolution, FileName)
% VelodromeModel creates a velodrome track black-line model that consists of 
% two straights, two circular arc bends and four transition curves between the 
% bends and the straights. The transition curves are based on two different
% Cesaro equations where the curvature along the arc length is defined. This 
% provides controlling of the centripetal acceleration during cornering and 
% has either G2 or G3 geometric continuity. 
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
%   path. This is a G3 continuous curve. 
% 
% The measurable features that define the track are:
%   L_L     The lap length. This is generally a known, fixed value. 
%   Y       The half-width between the two straights.
%   R       The turn radius at the bend apex.
%   If a clothoid is chosen
%       n   The power of the change in curvature with length. Usually 1 (linear)
% 
% The ratio of inputs Y/R has a limited range of feasible solutions dependent on
% L_L, R, and the curvature function that are checked before calculations begin. 
% 
% This function requires Matlab r2019b or later. 
% 
% Syntax
%   Track = VelodromeModel(Y, R)
%   Track = VelodromeModel(Y, R, 1)
%   Track = VelodromeModel(Y, R, 'sine')
%   Track = VelodromeModel(Y, R, n, L_L, Resolution, FileName)
% 
% Inputs
%   Y           (1 x 1 double)  [m] Half-width between the two straights.
%   R           (1 x 1 double)  [m] Radius of the circular bend arc. R < Y.
%   n       Two different options:
%               (1 x 1 double)  [-] Curvature exponent. n > 0. Default 1.
%                       OR
%               (1 x n char)    'sine' Use the sinusoidal curvature profile. 
%   L_L         (1 x 1 double)  [m] Lap length. Default 250.
%   Resolution  (1 x 1 double)  [m] Resolution of output data points. 
%                   Default 1. (Every 1 m along the datum line). 
%   FileName    (1 x m char)    Filename/path to save the output data.
%                   If FileName is not entered then the data will not be saved.
%                   File type based on extension: .txt, .dat, .csv, .xls, .xlsx
% 
% Output 
%   Data     	(m x 7 table) A table of the data vs track distance.
%                       m = L_L/Resolution + 1. 
%               Variables:
%                   Lap         [m]     Lap distance from the pursuit line
%                   X           [m]     x-coordinate
%                   Y           [m]     y-coordinate
%                   Curvature   [m^-1]  Curvature
%                   Radius      [m]     Radius of curvature
%                   dk_ds       [m^-2]  Cerivative of curvature w.r.t. Lap
%                   Tangent     [rad]   Tangential angle
%               Also has two structures 'Info' and 'Edge' in the table
%               custom properties that record calculation details. 
% 
% Example usage:
%   Track = VelodromeModel(23, 22, 1, 250, 0.25, 'TrackData.csv');
%   Track = VelodromeModel(23, 22, 'sine', 250, 0.25, 'TrackData.csv');
%   figure; plot(Track.X, Track.Y); axis equal
%   figure; plot(Track.Lap, Track.Curvature); 

%% Inputs 
arguments
    Y           (1,1) {double, mustBePositive}
    R           (1,1) {double, mustBePositive}
    n           (1,:) {double, char, string} = 1
    L_L         (1,1) {double, mustBePositive} = 250
    Resolution  (1,1) {double, mustBePositive} = 1
    FileName    (:,:) {char, string} = ''
end

nDataP = 2000; % [#] Number of data points for internal calculations

% Basic bounds
assert(Y < L_L/(2*pi), 'The half-width Y must be < L_L/(2*pi).')
assert(R < Y, 'The radius R must be < Y.')
assert(Resolution < L_L/25, 'The resolution must be << L_Lap.')

%% Transition curve calculations 
if isnumeric(n) 
    %% Curvature profile: Power equation 
    assert(isscalar(n) && n > 0, 'The exponent n must be scalar and > 0.')
    
    % Functions
    K  = @(v, L_T) v.^n/(R*L_T^n);              % Curvarture   
    Kd = @(v, L_T) n*v.^(n-1)/(R*L_T^n);        % Derivative 
    Ki = @(u, L_T) u.^(n+1)/(R*L_T^n*(n+1));    % Integral
    IC = @(t, L_T) integral(@(u) cos(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    IS = @(t, L_T) integral(@(u) sin(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    
    % Primary bounds for Y/R
    if R <= L_L/(2*pi*(n+1))
        A_Max = pi*R/2*(n+1);                   % From psi < pi/2
    else
        A_Max = (L_L - 2*pi*R)*(n+1)/(4*n);     % From L_S > 0
    end
    YonR_Max = K(A_Max, A_Max)*IS(A_Max, A_Max) + cos(Ki(A_Max, A_Max));
    assert(1 < Y/R && Y/R < YonR_Max, sprintf(...
        'Y/R (%.3f) must be in: 1 < Y/R < %.3f for R = %.3f & n = %g',...
        Y/R, YonR_Max, R, n))
    
    % Solving for A with the Newton–Raphson method.
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt(R*(Y - R)*(n+2)*(n+1));
    Er = 1;
    cc = 0;
    while abs(Er) > 1e-13
        % A1 = A0 - (IS(A0, A0) + R*cos(A0/R/(n+1)) - Y)/(sin(A0/R/(n+1))/(n+1));
        A1 = A0*(1-1/n) - (R*cos(A0/R/(n+1)) - Y)/(n/A0*IS(A0, A0));
        Er = A1 - A0;
        A0 = A1;
        cc = cc + 1;
        if cc > 2000, error('Calculation of A not converged'); end
    end
    A = A1;
    
    Style = sprintf('Power curvature, n: %g', n);
    Continuity = 'G2';
    
else
    %% Curvature profile: Sinusoidal 
    
    % Functions
    K  = @(v, L_T) (sin(pi/L_T*v - pi/2) + 1)/2/R;          % Curvarture 
    Kd = @(v, L_T)  cos(pi/L_T*v - pi/2)*pi/(2*R*L_T);      % Derivative 
    Ki = @(u, L_T) (-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R;  % Integral
    IC = @(t, L_T) integral(@(u) cos((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    IS = @(t, L_T) integral(@(u) sin((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    
    % Primary bounds for Y/R
    if R <= L_L/(4*pi)
        A_Max = pi*R;                               % From psi < pi/2
    else
        A_Max = (L_L - 2*pi*R)/2;                   % From L_S > 0
    end
    YonR_Max = K(A_Max, A_Max)*IS(A_Max, A_Max) + cos(Ki(A_Max, A_Max));
    assert(1 < Y/R && Y/R < YonR_Max, sprintf(...
        'Y/R (%.3f) must be in: 1 < Y/R < %.3f for R = %.3f',...
        Y/R, YonR_Max, R))
    
    % Solving for A with the Newton–Raphson method.
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt(4*pi*R*(Y - R)/(pi - 2));
    Er = 1;
    cc = 0;
    while abs(Er) > 1e-13
        A1 = A0 - (IS(A0, A0) + R*cos(A0/2/R) - Y)/(sin(A0/2/R)/2);
        Er = A1 - A0;
        A0 = A1;
        cc = cc + 1;
        if cc > 2000, error('Calculation of A not converged'); end
    end
    A = A1; 
    
    Style = 'Sinusoidal curvature';
    Continuity = 'G3';
    
end
L_T = A;

% Single-value results  
X       = IC(A, L_T) - R*sin(Ki(A, L_T));   % [m]   Bend centre X-coord
psi_1   = Ki(A, L_T);                       % [rad] Bend start tangent angle
theta   = pi/2 - psi_1;                     % [rad] Bend open angle
L_B     = theta*R;                          % [m]   Bend length
L_S     = L_L/4 - L_T - L_B;                % [m]   Straight length

% Per parametric distance t
t = linspace(0, A, nDataP)';                % [-]   Eqn parameter
for ii = 1:nDataP
    x(ii,1) = IC(t(ii), L_T);               % [m]   x coordinate
    y(ii,1) = IS(t(ii), L_T);               % [m]   y coordinate
end
s       = t;                                % [m]   Arc length
psi     = Ki(t, L_T);                       % [rad] Tangential angle
kappa   = K( t, L_T);                       % [m^-1] Curvature
dk_ds   = Kd(t, L_T);                       % [m^-2] Curvature derivative

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
XY.Str1 = [linspace(   0,  L_S, nDataP)', -Y*ones(nDataP, 1)];
XY.Str2 = [linspace( L_S, -L_S, nDataP)',  Y*ones(nDataP, 1)];
XY.Str3 = [linspace(-L_S,    0, nDataP)', -Y*ones(nDataP, 1)];

% Circular bends 
th = linspace(-theta, theta, nDataP)';         % [rad] Arc angle
XY.Bnd1 = [ xBc + R*cos(th),        yBc + R*sin(th)];
XY.Bnd2 = [-xBc + R*cos(th + pi),   yBc + R*sin(th + pi)];

%% Lap distance and subsequent parameters
%%%%% Lap distance 
Lap.Trn1 = L_S + s;
Lap.Trn2 = 1*L_S + 1*L_T + 2*L_B + s;
Lap.Trn3 = 3*L_S + 2*L_T + 2*L_B + s;
Lap.Trn4 = 3*L_S + 3*L_T + 4*L_B + s;

dsS = linspace(0, L_S, nDataP)';
Lap.Str1 = dsS;
Lap.Str2 = 1*L_S + 2*L_T + 2*L_B + 2*dsS;
Lap.Str3 = 3*L_S + 4*L_T + 4*L_B + 1*dsS;

dsC = R*linspace(0, 2*theta, nDataP)';
Lap.Bnd1 = 1*L_S + 1*L_T + dsC;
Lap.Bnd2 = 3*L_S + 3*L_T + 2*L_B + dsC;

%%%%% Curvature 
Curv.Trn1 = kappa;
Curv.Trn2 = flipud(kappa);
Curv.Trn3 = kappa;
Curv.Trn4 = flipud(kappa);

Curv.Str1 = zeros(nDataP, 1);
Curv.Str2 = zeros(nDataP, 1);
Curv.Str3 = zeros(nDataP, 1);

Curv.Bnd1 = ones(nDataP, 1)/R;
Curv.Bnd2 = ones(nDataP, 1)/R;

%%%%% Curvature derivative
CurvOnDs.Trn1 =  dk_ds;
CurvOnDs.Trn2 = -flipud(dk_ds);
CurvOnDs.Trn3 =  dk_ds;
CurvOnDs.Trn4 = -flipud(dk_ds);

CurvOnDs.Str1 = zeros(nDataP,1);
CurvOnDs.Str2 = zeros(nDataP,1);
CurvOnDs.Str3 = zeros(nDataP,1);

CurvOnDs.Bnd1 = zeros(nDataP,1);
CurvOnDs.Bnd2 = zeros(nDataP,1);

%%%%% Tangential angle 
Tangent.Trn1 = psi;
Tangent.Trn2 = -flipud(psi) + 2*theta + 2*psi_1;
Tangent.Trn3 = psi + pi;
Tangent.Trn4 = -flipud(psi) + 2*theta + 2*psi_1 + pi;

Tangent.Str1 = zeros(nDataP,1);
Tangent.Str2 = zeros(nDataP,1) +   pi;
Tangent.Str3 = zeros(nDataP,1) + 2*pi;

Tangent.Bnd1 = th + 1/2*pi;
Tangent.Bnd2 = th + 3/2*pi;

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

Comb.Curvature = [...
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
Info.n          = n;
Info.Y          = Y;
Info.X          = X;
Info.A          = A;
Info.R_Bnd      = R;
Info.L_Lap      = L_L;
Info.L_Str      = L_S;
Info.L_Trn      = L_T;
Info.L_Bnd      = L_B;
Info.theta      = theta;
Info.psi_1      = psi_1;
Info.BendCentre = [xBc, yBc];
Info.Trns       = cumsum([0, L_S, L_T, 2*L_B, L_T, 2*L_S, L_T, 2*L_B, L_T, L_S]);
Info.Resolution = Resolution;
Info.Continuity = Continuity;

%% Creating a consistently spaced table with the data
Track = table;
Track.Lap       = (0:Resolution:L_L)';
Track.X         = interp1(Comb.Lap, Comb.X,         Track.Lap, 'makima');
Track.Y         = interp1(Comb.Lap, Comb.Y,         Track.Lap, 'makima');
Track.Curvature = interp1(Comb.Lap, Comb.Curvature, Track.Lap, 'makima');
Track.Radius    = 1./Track.Curvature;
Track.dk_ds     = interp1(Comb.Lap, Comb.dk_ds,     Track.Lap, 'makima');
Track.Tangent   = interp1(Comb.Lap, Comb.Tangent,   Track.Lap, 'makima');

% Adjusting the tangential angle range to (-pi, +pi)
% This can be commented out if it is preferred to range in (0, 2*pi)
Track.Tangent(round(end/2):end) = Track.Tangent(round(end/2):end) - 2*pi;
% Floating-point error fix
Track.Tangent(abs(Track.Tangent)        < 1e-14) = 0;  
Track.Tangent(abs(Track.Tangent - pi)   < 1e-14) =  pi;
Track.Tangent(abs(Track.Tangent + pi)   < 1e-14) = -pi;

% Setting the track meta data
Track.Properties.VariableUnits = {'m', 'm', 'm', 'm^-1', 'm', 'm^-2', 'rad'};
Track = addprop(Track, {'Info', 'Edge'}, repmat({'table'},2,1)); 
Track.Properties.CustomProperties.Info = Info;
Track.Properties.CustomProperties.Edge = Edge;

% Saving the data to a file
if ~isempty(FileName)
    writetable(Track, FileName);
end
