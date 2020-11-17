% This script is to provide an example of using the VelodromeModel function

close all
clear
clc

% Inputs 
Y = 23.0;       % [m] Track half-width
R = 22.0;       % [m] Bend apex radius
L = 250;        % [m] Lap length
S = 0.1;        % [m] Resolution

% Power function
n = 1;          % [#] Curvature power

% Sinusoidal function
% n = 'sine';

Track = VelodromeModel(Y, R, n, L, S);

% To save the data 
% Track = VelodromeModel(Y, R, n, L, S, 'TrackData.csv');

% A structure with the track defining values 
Info = Track.Properties.CustomProperties.Info;

% A structure with the end coordinates of the different track segments 
Edge = Track.Properties.CustomProperties.Edge;

fprintf('%20s: %5.2f m\n', 'Straight length',   Info.L_Str)
fprintf('%20s: %5.2f m\n', 'Transition length', Info.L_Trn)
fprintf('%20s: %5.2f m\n', 'Bend length',   	Info.L_Bnd)
fprintf('Track:\n')
disp(head(Track))

figure; 

%%%%%%%%%% (x, y) coordinates
subplot(3,2,[1,2])
hold on
box  on
plot(Track.X, Track.Y)
scatter(Edge.Bnd_Centre(:,1), Edge.Bnd_Centre(:,2), 'ok')
scatter(Edge.Str_Edge(:,1),   Edge.Str_Edge(:,2),   'sk')
scatter(Edge.Bnd_Edge(:,1),   Edge.Bnd_Edge(:,2),   'sk')
axis equal
axis tight
xlabel('X [m]')
ylabel('Y [m]')
title('Plan view')
set(gca, 'XLim',get(gca, 'XLim')*1.05)
set(gca, 'YLim',get(gca, 'YLim')*1.10)

%%%%%%%%%% Curvature
subplot(3,2,3)
plot(Track.Lap, Track.Curvature)
xlabel('Lap Position [m]')
ylabel('\kappa [m^{-1}]')
title('Curvature')
xlim([0, L])

%%%%%%%%%% Turn radius
subplot(3,2,4)
plot(Track.Lap, Track.Radius)
xlabel('Lap Position [m]')
ylabel('\rho [m]')
title('Turn radius')
xlim([0, L])
ylim([0, 100])

%%%%%%%%%% Derivative of curvature
subplot(3,2,5)
plot(Track.Lap, Track.dk_ds)
xlabel('Lap Position [m]')
ylabel('d\kappa/ds [m^{-2}]')
title('Derivative of curvature')
xlim([0, L])

%%%%%%%%%% Tangent
subplot(3,2,6)
plot(Track.Lap, Track.Tangent)
xlabel('Lap Position [m]')
ylabel('\psi [rad]')
title('Tangential angle')
xlim([0, L])
