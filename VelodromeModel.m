function Track = VelodromeModel(Y, R, n, L_L, Opts)
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
%   Track = VelodromeModel(Y, R)
%   Track = VelodromeModel(Y, R, 1)
%   Track = VelodromeModel(Y, R, 'sine')
%   Track = VelodromeModel(Y, R, n, L_L)
%   Track = VelodromeModel(..., 'Bank',[<value>, <value>])
%   Track = VelodromeModel(..., 'Width',<value>)
%   Track = VelodromeModel(..., 'Resolution',<value>)
%   Track = VelodromeModel(..., 'FileName',<name>)
% 
% Inputs (Ordered)
%   Y           (1 x 1 double)  [m] Half-span between the two straights.
%   R           (1 x 1 double)  [m] Radius of the circular bend arc. R < Y.
%   n       Two different options:
%               (1 x 1 double)  [-] Curvature exponent. n > 0. Default 1.
%                       OR
%               (1 x n char)    'sine' Use the sinusoidal curvature profile. 
%   L_L         (1 x 1 double)  [m] Lap length. Default 250.
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
% Output 
%   Track     	(m x 9+ table) A table of the data vs track distance.
%                       m = L_L/Resolution + 1. 
%               Variables:
%                   Lap         [m]     Lap distance from the pursuit line
%                   X           [m]     x-coordinate of the datum line
%                   Y           [m]     y-coordinate of the datum line
%                   Z           [m]     z-coordinate of the datum line (0)
%                   Curvature   [m^-1]  Curvature
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
% 
% Example usage:
%   Track = VelodromeModel(23, 22);
%   Track = VelodromeModel(23, 22, 'sine', 250, 'Bank',[13, 43], 'Width',7.5);
%   Track = VelodromeModel(23, 22,    1,   250, 'Resolution',0.25);
%   Track = VelodromeModel(23, 22,    1,   250, 'FileName','TrackData.csv');
%   figure; stackedplot(Track, 'XVariable','Lap');
%   figure; plot(Track.X, Track.Y); axis equal;
%   figure; plot(Track.Lap, Track.Curvature); 
%
% Submitted to Sports Engineering 
% 'The impact of transition curve design on the accuracy of velodrome models' 

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
if isnumeric(n) 
    %% Curvature profile: Power equation (17) 
    assert(isscalar(n) && n > 0, 'The exponent n must be scalar and > 0.')
    
    Style = sprintf('Power, n: %g', n);
    Continuity = 'G2';
    
    % Functions
    K  = @(v, L_T) v.^n/(R*L_T^n);                              % Curvature   
    Kd = @(v, L_T) n*v.^(n-1)/(R*L_T^n);                        % 1st derivative 
    K2d= @(v, L_T) n*(n-1)*v.^(n-2)/(R*L_T^n);                  % 2nd derivative 
    Ki = @(u, L_T) u.^(n+1)/(R*L_T^n*(n+1));                    % Integral
    IC = @(t, L_T) integral(@(u) cos(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    IS = @(t, L_T) integral(@(u) sin(u.^(n+1)/(R*L_T^n*(n+1))), 0, t);
    
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt((n+2)*(n+1)*R*(Y - R));
    
elseif strcmpi(n(1), 's')
    %% Curvature profile: Sinusoidal (18) 
    Style = 'Sinusoidal';
    Continuity = 'G3';
    
    % Functions
    K  = @(v, L_T) (sin(pi/L_T*v - pi/2) + 1)/2/R;              % Curvature 
    Kd = @(v, L_T)  cos(pi/L_T*v - pi/2)*pi/(2*R*L_T);          % 1st derivative 
    K2d= @(v, L_T) -sin(pi/L_T*v - pi/2)*pi^2/(2*R*L_T^2);      % 2nd derivative 
    Ki = @(u, L_T) (-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R;      % Integral
    IC = @(t, L_T) integral(@(u) cos((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    IS = @(t, L_T) integral(@(u) sin((-L_T/pi*cos(pi/L_T*u - pi/2) + u)/2/R), 0, t);
    n  = 1;
    
    % The initial value is found with the first-degree Maclaurin polynomial
    A0 = sqrt(4*pi/(pi - 2)*R*(Y - R));
    
else
    error('Unknown model type: ''%s''.', n)
    
end

% Primary bounds for Y/R
if R <= L_L/(2*pi*(n+1))
    A_Max = pi*R/2*(n+1);                   % [-]    (15) (psi < pi/2)
else
    A_Max = (L_L - 2*pi*R)*(n+1)/(4*n);     % [-]    (16) (L_S > 0)
end
YonR_Max = K(A_Max, A_Max)*IS(A_Max, A_Max) + cos(Ki(A_Max, A_Max));
assert(1 < Y/R && Y/R < YonR_Max, sprintf(...
    'Y/R (%.3f) must be in: (1 < Y/R < %.4f) for ''%s'' curvature with R = %.4g.',...
    Y/R, YonR_Max, Style, R))

% Solving for A with the Newton–Raphson method.
Er = 1;
cc = 0;
while abs(Er) > 1e-13
    A1 = A0 - (IS(A0, A0) + R*cos(A0/R/(n+1)) - Y)/(sin(A0/R/(n+1))/(n+1));
    Er = A1 - A0;
    A0 = A1;
    cc = cc + 1;
    if cc > 2000, error('Calculation of A not converged'); end
end
A   = A1;                                   % [m]    ( 9) tau_1 
L_T = A1;                                   % [m]    (10) Transition length

% Single-value results  
X       = IC(A, L_T) - R*sin(Ki(A, L_T));   % [m]    ( 9) Bend centre X-coord
psi_1   = Ki(A, L_T);                       % [rad]  ( 3) Bend start tangent 
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
t = A*ArrayA;                               % [-]    Eqn parameter
for ii = 1:nDataP
    x(ii,1) = IC(t(ii), L_T);               % [m]    ( 1) x coordinate
    y(ii,1) = IS(t(ii), L_T);               % [m]    ( 1) y coordinate
end
s       = t;                                % [m]    ( 2) Arc length
psi     = Ki(t, L_T);                       % [rad]  ( 3) Tangential angle
kappa   = K( t, L_T);                       % [m^-1] ( 4) Curvature
dk_ds   = Kd(t, L_T);                       % [m^-2] ( 5) Curvature 1st derivative
d2k_ds2 = K2d(t, L_T);                      % [m^-3]      Curvature 2nd derivative
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
Lap.Trn1 = L_S + s;
Lap.Trn2 = 1*L_S + 1*L_T + 2*L_B + s;
Lap.Trn3 = 3*L_S + 2*L_T + 2*L_B + s;
Lap.Trn4 = 3*L_S + 3*L_T + 4*L_B + s;

Lap.Str1 = s_S;
Lap.Str2 = 1*L_S + 2*L_T + 2*L_B + 2*s_S;
Lap.Str3 = 3*L_S + 4*L_T + 4*L_B + 1*s_S;

Lap.Bnd1 = 1*L_S + 1*L_T + R*th;
Lap.Bnd2 = 3*L_S + 3*L_T + 2*L_B + R*th;

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

%%%%% Curvature 1st derivative
CurvOnDs.Trn1 =  dk_ds;
CurvOnDs.Trn2 = -flipud(dk_ds);
CurvOnDs.Trn3 =  dk_ds;
CurvOnDs.Trn4 = -flipud(dk_ds);

CurvOnDs.Str1 = zeros(nDataP,1);
CurvOnDs.Str2 = zeros(nDataP,1);
CurvOnDs.Str3 = zeros(nDataP,1);

CurvOnDs.Bnd1 = zeros(nDataP,1);
CurvOnDs.Bnd2 = zeros(nDataP,1);

%%%%% Curvature 2nd derivative
CurvOnDs2.Trn1 =  d2k_ds2;
CurvOnDs2.Trn2 = -flipud(d2k_ds2);
CurvOnDs2.Trn3 =  d2k_ds2;
CurvOnDs2.Trn4 = -flipud(d2k_ds2);

%%%%% Tangential angle 
Tangent.Trn1 = psi;
Tangent.Trn2 = -flipud(psi) + 2*theta + 2*psi_1;
Tangent.Trn3 = psi + pi;
Tangent.Trn4 = -flipud(psi) + 2*theta + 2*psi_1 + pi;

Tangent.Str1 = zeros(nDataP,1);
Tangent.Str2 = zeros(nDataP,1) +   pi;
Tangent.Str3 = zeros(nDataP,1) + 2*pi;

Tangent.Bnd1 = th - theta + 1/2*pi;
Tangent.Bnd2 = th - theta + 3/2*pi;

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
Info.Transition = Transition;
Info.Continuity = Continuity;
Info.Bank       = Opts.Bank;
Info.Width      = Opts.Width;
Info.Resolution = Opts.Resolution;

%% Creating a consistently spaced table with the data
Track = table;
Track.Lap       = (0:Opts.Resolution:L_L)';
Track.X         = interp1(Comb.Lap, Comb.X,         Track.Lap, 'makima');
Track.Y         = interp1(Comb.Lap, Comb.Y,         Track.Lap, 'makima');
Track.Z         = zeros(height(Track),1);
Track.Curvature = interp1(Comb.Lap, Comb.Curvature, Track.Lap, 'makima');
Track.Radius    = 1./Track.Curvature;
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
Track.Properties.VariableUnits = {'m','m','m','m', 'm^-1','m','m^-2','m^-3','rad'};
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
