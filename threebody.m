%Lecture 18

%three body sim courtesy of dr.kamrin
clear all
close all
clc

% I have chosen to normalize everything by G, so I will set G=1.
global m1 m2 m3 G
G = 1;

%% --- Mass Profiles (Uncomment the case you want to run) ---

% Case 1: One-body limit, m2 >> m1. Keplerian orbits
% m1 = 1;
% m2 = 100;
% m3 = 0;

% Case 2: Two-body problem, m1 similar to m2. Stable orbits
% m1 = 9;
% m2 = 7;
% m3 = 0;

% Case 3: Three-body problem. Semi-stable epoch
% m1 = 7; 
% m2 = 6;
% m3 = 8;

% Case 4: Three-body problem. Chaos
m1 = 7; 
m2 = 8;
m3 = 3;

%% --- Initial Conditions and ODE Solver ---
% Initial state vector format: [v1; v2; v3; r1; r2; r3]
initial_state = [3/m1; 0; 0; ...   % v1 (x, y, z)
                 -3/m2; 0; 0; ...  % v2 (x, y, z)
                 0; 0; 0; ...      % v3 (x, y, z)
                 1; 1; 0; ...      % r1 (x, y, z)
                 2; 3; 4; ...      % r2 (x, y, z)
                 5; 1; 2];         % r3 (x, y, z)

tspan = 0:0.01:10;

% Solve the differential equations using ode15s
[t, Uout] = ode15s(@threebody, tspan, initial_state);

% Extract positions and velocities from the solver output matrix
v1 = Uout(:, 1:3);
v2 = Uout(:, 4:6);
v3 = Uout(:, 7:9);
r1 = Uout(:, 10:12);
r2 = Uout(:, 13:15);
r3 = Uout(:, 16:18);

%% --- Dynamically Set Up Plot Axes Limits ---
if m3 ~= 0
    xmin = min([r1(:,1); r2(:,1); r3(:,1)]);
    xmax = max([r1(:,1); r2(:,1); r3(:,1)]);
    ymin = min([r1(:,2); r2(:,2); r3(:,2)]);
    ymax = max([r1(:,2); r2(:,2); r3(:,2)]);
    zmin = min([r1(:,3); r2(:,3); r3(:,3)]);
    zmax = max([r1(:,3); r2(:,3); r3(:,3)]);
else
    xmin = min([r1(:,1); r2(:,1)]);
    xmax = max([r1(:,1); r2(:,1)]);
    ymin = min([r1(:,2); r2(:,2)]);
    ymax = max([r1(:,2); r2(:,2)]);
    zmin = min([r1(:,3); r2(:,3)]);
    zmax = max([r1(:,3); r2(:,3)]);
end

%% --- Animation Loop ---
figure('Name', 'Three-Body Problem Simulation', 'Color', 'w');

for i = 1:length(t)
    % Body 1: Current position and historical trajectory trace
    plot3(r1(i,1), r1(i,2), r1(i,3), '*', 'MarkerSize', 20)
    hold on
    plot3(r1(1:i,1), r1(1:i,2), r1(1:i,3), 'b-')
    
    % Body 2: Current position and historical trajectory trace
    plot3(r2(i,1), r2(i,2), r2(i,3), '*', 'MarkerSize', 20)
    plot3(r2(1:i,1), r2(1:i,2), r2(1:i,3), 'r-')
    
    % Body 3: Only process and plot if mass 3 is present
    if m3 ~= 0
        plot3(r3(i,1), r3(i,2), r3(i,3), '*', 'MarkerSize', 20)
        plot3(r3(1:i,1), r3(1:i,2), r3(1:i,3), 'g-')
    end
    
    % Formatting the plot environment
    axis([xmin, xmax, ymin, ymax, zmin, zmax])
    grid on
    view(3) % Forces a 3D perspective angle
    xlabel('X Space Axis')
    ylabel('Y Space Axis')
    zlabel('Z Space Axis')
    title(sprintf('Time: %.2f s', t(i)))
    
    pause(0.0001)
    hold off
end

%% derivative calculations
function Udot = threebody(~, U)
    global m1 m2 m3 G

    % Extract individual 3D vectors out of the long vector U
    v1 = U(1:3);
    v2 = U(4:6);
    v3 = U(7:9);
    r1 = U(10:12);
    r2 = U(13:15);
    r3 = U(16:18);

    % The time-rate of change of each velocity vector is given by Newton's Law of Gravity
    v1dot = G * ((m2 / norm(r2 - r1)^3) * (r2 - r1) + (m3 / norm(r3 - r1)^3) * (r3 - r1));
    v2dot = G * ((m1 / norm(r1 - r2)^3) * (r1 - r2) + (m3 / norm(r3 - r2)^3) * (r3 - r2));
    v3dot = G * ((m1 / norm(r1 - r3)^3) * (r1 - r3) + (m2 / norm(r2 - r3)^3) * (r2 - r3));

    % The time-rate of change of each position vector is simply its velocity
    r1dot = v1;
    r2dot = v2;
    r3dot = v3;

    % Reassemble the derivatives back into a single column vector Udot
    Udot = [v1dot; v2dot; v3dot; r1dot; r2dot; r3dot];
end
