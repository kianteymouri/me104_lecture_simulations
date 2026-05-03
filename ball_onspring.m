%%Lecture 8 example: code courtesy of Dr.Kamrin

clear; clc; close all;
global k m l0

%% 1. Parameters & Initial Conditions
vx0 = 1; vy0 = 0.5; x0 = 1.5; % Added vy0 to ensure orbital motion
k = 100; m = 1; l0 = 1.8;
Tmax = 5; Tinc = 1/300;

%% 2. The ODE Solver (This creates 'Usol')
% This must happen BEFORE you reference 'Usol' in the next lines
[t, Usol] = ode45(@springexample, [0:Tinc:Tmax], [x0, vx0, 0, vy0/x0]);

%% 3. Visual Setup (Now 'Usol' exists!)
max_r = max(Usol(:,1)) * 1.3; 
figure('Color', 'w');

%% 4. Animation Loop
for i = 1:10:length(t)
    % Convert polar results back to Cartesian for plotting
    x = Usol(i,1) * cos(Usol(i,3));
    y = Usol(i,1) * sin(Usol(i,3));
    
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 
    hold on;
    plot([0, x], [0, y], 'k-', 'LineWidth', 2); % Spring
    
    % Fixed the LineWidth syntax from your previous version
    plot(Usol(1:i,1).*cos(Usol(1:i,3)), Usol(1:i,1).*sin(Usol(1:i,3)), ...
        'LineStyle', ':', 'LineWidth', 1.5, 'Color', 'b'); 
        
    axis([-max_r max_r -max_r max_r]);
    axis square; grid on;
    hold off;
    drawnow;
end

%% 5. ODE Function Definition
function Udot = springexample(t, U)
    global k m l0
    r = U(1); alpha = U(2); theta = U(3); beta = U(4);
    
    % Equations of motion derived from F = ma
    rdot = alpha;
    alphadot = -(k/m)*(r-l0) + r*beta^2; % Radial acceleration
    thetadot = beta;
    betadot = -2*alpha*beta/r;           % Tangential acceleration
    
    Udot = [rdot; alphadot; thetadot; betadot];
end
