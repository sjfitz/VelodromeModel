% This script is to provide an example of using the VelodromeModel function

close all
clear
clc

% Inputs 
Y = 23.0;       % [m]   Track half-span
R = 22.0;       % [m]   Bend apex radius
L = 250;        % [m]   Lap length
S = 0.1;        % [m]   Resolution
Bank_min = 13;  % [deg] Minimum bank angle
Bank_max = 43;  % [deg] Maximum bank angle
Width = 7.5;    % [m]   Track width

Bank = [Bank_min, Bank_max];

% Power function
% n = 1;          % [#]   Curvature power

% Sinusoidal function
n = 'sine';

% Simple inputs 
Track = VelodromeModel(Y, R, n);

% All inputs
% Track = VelodromeModel(Y, R, n, L, 'Bank',Bank, 'Width',Width, 'Resolution',S);

% To save the data 
% Track = VelodromeModel(Y, R, n, L, 'FileName','TrackData.csv');

% A structure with the track defining values 
Info = Track.Properties.CustomProperties.Info;

% A structure with the end coordinates of the different track segments 
Edge = Track.Properties.CustomProperties.Edge;

fprintf('%20s: %5.2f m\n', 'Straight length',   Info.L_Str)
fprintf('%20s: %5.2f m\n', 'Transition length', Info.L_Trn)
fprintf('%20s: %5.2f m\n', 'Bend length',   	Info.L_Bnd)
fprintf('Track:\n')
disp(head(Track))

%% Plotting - All data 
figure; 

%%%%%%%%%% (x, y) coordinates
subplot(3,2,1)
hold on
box  on
plot(Track.X, Track.Y)
scatter(Edge.Bnd_Centre(:,1), Edge.Bnd_Centre(:,2), 'ok')
scatter(Edge.Str_Edge(:,1),   Edge.Str_Edge(:,2),   'sk')
scatter(Edge.Bnd_Edge(:,1),   Edge.Bnd_Edge(:,2),   'sk')
axis equal
xlabel('X [m]')
ylabel('Y [m]')
title('Plan view')
axis tight
set(gca, 'XLim',get(gca, 'XLim')*1.05)
set(gca, 'YLim',get(gca, 'YLim')*1.10)
% axis padded

%%%%%%%%%% Bank Angle
if ~isequal(Info.Bank, [0, 0])
    subplot(3,2,2)
    plot(Track.Lap, Track.BankAngle)
    xlabel('Lap Position [m]')
    ylabel('\beta [deg]')
    title('Bank Angle')
    xlim([0, L])
end

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
ylabel('\kappa'' [m^{-2}]')
title('Derivative of curvature')
xlim([0, L])

%%%%%%%%%% Tangent
subplot(3,2,6)
plot(Track.Lap, Track.Tangent)
xlabel('Lap Position [m]')
ylabel('\psi [rad]')
title('Tangential angle')
xlim([0, L])
yticks((-1:0.5:1)*pi)
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})

%% Plotting - 3D view 
if Info.Width ~= 0
    figure;
    hold on
    box  on
    grid on
    
    plot3(Track.X, Track.Y, Track.Z, 'k', 'linewidth',2);
    plot3(Track.X_Top, Track.Y_Top, Track.Z_Top, 'k', 'linewidth',2);
    surf(...
        [Track.X, Track.X_Top], ...
        [Track.Y, Track.Y_Top], ...
        [Track.Z, Track.Z_Top], 'EdgeColor','none');
    
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title('3D view')
    view([15, 20])
end
