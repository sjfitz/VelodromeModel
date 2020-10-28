function Track = VelodromeModel(Y, R, n, L_L, Resolution, FileName)
% VelodromeModel creates a velodrome track black-line model that consists of 
% two straights, two circular arc bends and four transition curves between the 
% bends and the straights. The transition curves follow the form of a clothoid
% (also known as a Euler spiral or Cornu spiral). A clothoid has the unique 
% property that its curvature is a polynomial (usually linear) function of its 
% arc length. This provides controlling of the centripetal acceleration during 
% cornering and has G2 geometric continuity. 
% 
% The features that define the track are:
%   L_L     The lap length. This is generally a known, fixed value. 
%   Y       The half-width between the two straights.
%   R       The turn radius at the bend apex.
%   n       The power of the change in curvature with length. Usually 1 (linear)
% The ratio of inputs Y/R has a limited range of feasible solutions dependent on
% L_L, R, & n that are checked before calculations begin. 
% 
% Track = VelodromeModel(Y, R, n, L_L, Resolution, FileName)
% 
% Example usage:
%   Track = VelodromeModel(22, 21, 1, 250, 0.25, 'TrackData.csv')
%   figure; plot(Track.X, Track.Y); axis equal
%   figure; plot(Track.Lap, Track.Curvature); 
% 
% Inputs
%   Y           (1 x 1 double) [m]      Half-width between the two straights.
%                   Limits: Y < L_L/(2*pi). See paper for further limits. 
%   R           (1 x 1 double) [m]      Radius of the circular bend arc.
%                   Limits: R < Y.          See paper for further limits. 
%   n           (1 x 1 double) [-]      Curvature exponent. 
%                   Limits: n > 0.          Default 1.
%   L_L         (1 x 1 double) [m]      Lap length. 
%                                           Default 250.
%   Resolution  (1 x 1 double) [m]      Resolution of output data points. 
%                   Default 1. (Every 1 m along the datum line). 
%   FileName    (1 x m char)            Filename/path to save the output data.
%                   If FileName is not entered then the data will not be saved.
%                   File type based on extension: .txt, .dat, .csv, .xls, .xlsx
% 
% Output 
%   Data     	(m x 7 table) A table of the data vs track distance.
%                       m = L_L/Resolution + 1. 
%               Variables:
%                   Lap         [m]     Lap distance from the pursuit line
%                   Curvature   [m^-1]  Track curvature
%                   Radius      [m]     Track radius
%                   X           [m]     x-coordinate
%                   Y           [m]     y-coordinate
%                   dkOnds      [m^-2]  Track curvature derivative w.r.t. Lap
%                   Tangent     [rad]   Tangential angle
%               Also has two structures 'Info' and 'Edge' in the table
%               custom properties that records calculation details. 
% 
% Shaun Fitzgerald 
% Created 2020-10-25 

%% Inputs 
arguments
    Y           (1,1) double
    R           (1,1) double
    n           (1,1) double = 1
    L_L         (1,1) double = 250
    Resolution  (1,1) double = 1
    FileName    (:,:) {string, char} = ''
end

% Basic bounds
assert(R < Y, 'The radius R must be < Y')
assert(n > 0, 'The exponent n must be > 0')
assert(Resolution < L_L/25, 'The resolution must be much less than the lap')
assert(2*pi*Y < L_L, 'The half-width Y must be < L_L/(2*pi).')

nDataP = 1000;     % [#] Number of data points for internal calculations

%% Clothoid calculations 
% Primary bounds for Y/R - Checking if the solution is feasible  
if R <= L_L/(2*pi*(n+1))    % From psi < pi/2
    A_Max = (pi*(n+1)/2)^(1/(n+1));
else                        % From L_S > 0
    A_Max = ((L_L/4/R - pi/2)*(n+1)/n)^(1/(n+1));
end
YonR_Max = A_Max^n*IS(A_Max, n) + cos(A_Max^(n+1)/(n+1));
assert(1 < Y/R && Y/R < YonR_Max, sprintf(...
    'Y/R (%.3f) must be in the range: 1 < Y/R < %.3f for n = %g & R = %.3f',...
    Y/R, YonR_Max, n, R))

% Solving for A with the Newton–Raphson method.
% The starting value is found by solving the first-degree Maclaurin polynomial
A0 = ((Y/R-1)*(n+2)*(n+1))^(1/2/(n+1)); 
Er = 1;
cc = 0;
while abs(Er) > 1e-13
    A1 = A0*(1-1/n) - (R*cos(A0^(n+1)/(n+1)) - Y)/(n*R*A0^(n-1)*IS(A0, n));
    Er = A1 - A0;
    A0 = A1;
    cc = cc + 1;
    if cc > 1000, error('Calculation of A not converged'); end
end
A = A1;

% Other calculations 
a   = R*A^n;                                    % [m]   Scale factor
X   = R*A^n*IC(A, n) - R*sin(A^(n+1)/(n+1));    % [m]   X-coord of the end bend
psi = A^(n+1)/(n+1);                            % [rad] Tangential angle 
theta = pi/2 - psi;                             % [rad] Bend open angle
L_T = R*A^(n+1);                                % [m]   Transition length
L_B = theta*R;                                  % [m]   Bend length 
L_S = L_L/4 - L_T - L_B;                        % [m]   Straight length 

% Per parametric distance t
t = linspace(0, A, nDataP)';                    % [-]   Eqn parameter
% This can be optionally replaced with a parfor loop for a speed improvement.
for ii = 1:nDataP
    x(ii,1) = a*IC(t(ii), n);                   % [m]   x coordinate
    y(ii,1) = a*IS(t(ii), n);                   % [m]   y coordinate
end
s       = R*A^n*t;                              % [m]    Arc length
kappa   = t.^n/R/A^n;                           % [m^-1] Curvature
dkOnds  = n*t.^(n-1)/R^2/A^(2*n);               % [m^-2] Curvature derivative
psi_t   = t.^(n+1)/(n+1);                       % [rad]  Tangential angle

% Moving the origin to the centre of the velodrome
xT  = x + L_S;
yT  = y - Y;
xBc = X + L_S;
yBc = 0;

%% The (x, y) coordinates for each track segment
% Edge points of the straights 
Edge.Str_Edge(1,:) = [ L_S, -Y];
Edge.Str_Edge(2,:) = [ L_S,  Y]; 
Edge.Str_Edge(3,:) = [-L_S,  Y]; 
Edge.Str_Edge(4,:) = [-L_S, -Y]; 

% Edge points of the circular bend arcs
Edge.Cba_Edge(1,:) = [ xBc + R*cos(theta), yBc - R*sin(theta)];
Edge.Cba_Edge(2,:) = [ xBc + R*cos(theta), yBc + R*sin(theta)];
Edge.Cba_Edge(3,:) = [-xBc - R*cos(theta), yBc + R*sin(theta)];
Edge.Cba_Edge(4,:) = [-xBc - R*cos(theta), yBc - R*sin(theta)];

% Centre of the circular bend arc
Edge.Cba_Centre(1,:) = [ xBc, yBc];
Edge.Cba_Centre(2,:) = [-xBc, yBc];

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
th   = linspace(-theta, theta, nDataP)';         % [rad] Arc angle
XY.Cba1 = [ xBc + R*cos(th),        yBc + R*sin(th)];
XY.Cba2 = [-xBc + R*cos(th + pi),   yBc + R*sin(th + pi)];

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
Lap.Cba1 = 1*L_S + 1*L_T + dsC;
Lap.Cba2 = 3*L_S + 3*L_T + 2*L_B + dsC;

%%%%% Curvature 
Curv.Trn1 = kappa;
Curv.Trn2 = flipud(kappa);
Curv.Trn3 = kappa;
Curv.Trn4 = flipud(kappa);

Curv.Str1 = zeros(nDataP, 1);
Curv.Str2 = zeros(nDataP, 1);
Curv.Str3 = zeros(nDataP, 1);

Curv.Cba1 = ones(nDataP, 1)/R;
Curv.Cba2 = ones(nDataP, 1)/R;

%%%%% Curvature derivative
CurvOnDs.Trn1 =  dkOnds;
CurvOnDs.Trn2 = -flipud(dkOnds);
CurvOnDs.Trn3 =  dkOnds;
CurvOnDs.Trn4 = -flipud(dkOnds);

CurvOnDs.Str1 = zeros(nDataP,1);
CurvOnDs.Str2 = zeros(nDataP,1);
CurvOnDs.Str3 = zeros(nDataP,1);

CurvOnDs.Cba1 = zeros(nDataP,1);
CurvOnDs.Cba2 = zeros(nDataP,1);

%%%%% Tangential angle 
Tangent.Trn1 = psi_t;
Tangent.Trn2 = -flipud(psi_t) + 2*theta + 2*psi;
Tangent.Trn3 = psi_t + pi;
Tangent.Trn4 = -flipud(psi_t) + 2*theta + 2*psi + pi;

Tangent.Str1 = zeros(nDataP,1);
Tangent.Str2 = zeros(nDataP,1) +   pi;
Tangent.Str3 = zeros(nDataP,1) + 2*pi;

Tangent.Cba1 = th + 1/2*pi;
Tangent.Cba2 = th + 3/2*pi;

%% Combining the data into one (unequally spaced) lap distance 
Comb.Lap = [...
    Lap.Str1(1:end-1); ...
    Lap.Trn1(1:end-1); ...
    Lap.Cba1(1:end-1); ...
    Lap.Trn2(1:end-1); ...
    Lap.Str2(1:end-1); ...
    Lap.Trn3(1:end-1); ...
    Lap.Cba2(1:end-1); ...
    Lap.Trn4(1:end-1); ...
    Lap.Str3(1:end-1)];

Comb.Curvature = [...
    Curv.Str1(1:end-1); ...
    Curv.Trn1(1:end-1); ...
    Curv.Cba1(1:end-1); ...
    Curv.Trn2(1:end-1); ...
    Curv.Str2(1:end-1); ...
    Curv.Trn3(1:end-1); ...
    Curv.Cba2(1:end-1); ...
    Curv.Trn4(1:end-1); ...
    Curv.Str3(1:end-1)];

Comb.X = [...
    XY.Str1(1:end-1,1); ...
    XY.Trn1(1:end-1,1); ...
    XY.Cba1(1:end-1,1); ...
    XY.Trn2(1:end-1,1); ...
    XY.Str2(1:end-1,1); ...
    XY.Trn3(1:end-1,1); ...
    XY.Cba2(1:end-1,1); ...
    XY.Trn4(1:end-1,1); ...
    XY.Str3(1:end-1,1)];

Comb.Y = [...
    XY.Str1(1:end-1,2); ...
    XY.Trn1(1:end-1,2); ...
    XY.Cba1(1:end-1,2); ...
    XY.Trn2(1:end-1,2); ...
    XY.Str2(1:end-1,2); ...
    XY.Trn3(1:end-1,2); ...
    XY.Cba2(1:end-1,2); ...
    XY.Trn4(1:end-1,2); ...
    XY.Str3(1:end-1,2)];

Comb.dkOnds = [...
    CurvOnDs.Str1(1:end-1); ...
    CurvOnDs.Trn1(1:end-1); ...
    CurvOnDs.Cba1(1:end-1); ...
    CurvOnDs.Trn2(1:end-1); ...
    CurvOnDs.Str2(1:end-1); ...
    CurvOnDs.Trn3(1:end-1); ...
    CurvOnDs.Cba2(1:end-1); ...
    CurvOnDs.Trn4(1:end-1); ...
    CurvOnDs.Str3(1:end-1)];

Comb.Tangent = [...
    Tangent.Str1(1:end-1); ...
    Tangent.Trn1(1:end-1); ...
    Tangent.Cba1(1:end-1); ...
    Tangent.Trn2(1:end-1); ...
    Tangent.Str2(1:end-1); ...
    Tangent.Trn3(1:end-1); ...
    Tangent.Cba2(1:end-1); ...
    Tangent.Trn4(1:end-1); ...
    Tangent.Str3(1:end-1)];

%% Saving the primary information
Info.Resolution = Resolution;
Info.n          = n;
Info.Y          = Y;
Info.X          = X;
Info.R_cba      = R;
Info.L_Lap      = L_L;
Info.L_Str      = L_S;
Info.L_Trn      = L_T;
Info.L_cba      = L_B;
Info.theta      = theta;
Info.psi        = psi;
Info.a          = a;
Info.A          = A;
Info.BendCentre = [xBc, yBc];

%% Creating a consistently spaced table with the data
Track = table;
Track.Lap       = (0:Resolution:L_L)';
Track.X         = interp1(Comb.Lap, Comb.X,         Track.Lap, 'makima');
Track.Y         = interp1(Comb.Lap, Comb.Y,         Track.Lap, 'makima');
Track.Curvature = interp1(Comb.Lap, Comb.Curvature, Track.Lap, 'makima');
Track.Radius    = 1./Track.Curvature;
Track.dkOnds    = interp1(Comb.Lap, Comb.dkOnds,    Track.Lap, 'makima');
Track.Tangent   = interp1(Comb.Lap, Comb.Tangent,   Track.Lap, 'makima');

% Adjusting the tangential angle range to (-pi, +pi)
% This can be commented out if it is preferred to range in (0, 2*pi)
Track.Tangent(round(end/2):end) = Track.Tangent(round(end/2):end) - 2*pi;
Track.Tangent(abs(Track.Tangent) < 1e-14) = 0;  % Floating-point error fix
Track.Tangent(abs(Track.Tangent - pi) < 1e-14) =  pi;
Track.Tangent(abs(Track.Tangent + pi) < 1e-14) = -pi;

% Setting the track meta data
Track.Properties.VariableUnits = {'m', 'm', 'm', 'm^-1', 'm', 'm^-2', 'rad'};
Track = addprop(Track, {'Info', 'Edge'}, repmat({'table'},2,1)); 
Track.Properties.CustomProperties.Info = Info;
Track.Properties.CustomProperties.Edge = Edge;

% Saving the data to a file
if ~isempty(FileName)
    writetable(Track, FileName);
end

end

% Clothoid integrals. These could be anonymous functions, but this is ~4% faster
function Int = IC(t, n), Int = integral(@(x) cos(x.^(n+1)/(n+1)), 0, t); end
function Int = IS(t, n), Int = integral(@(x) sin(x.^(n+1)/(n+1)), 0, t); end
