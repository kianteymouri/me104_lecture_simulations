% ME 104 Lecture 2: Car being pulled up a wedge
clear; clc; close all;

%% 1. Parameters
v_rope = 1.0;       %arb speed being pulled up
phi_deg = 30;       %rope angle
theta_deg = 20;     %wedge angle
L_total = 12;       %arb length of rope

phi = deg2rad(phi_deg);
theta = deg2rad(theta_deg);

t_max = 5;          
fps = 30;
t = linspace(0, t_max, t_max * fps);

%%equations from lecture
% l1 + l2 = constant
l2_0 = 1;
l1_0 = L_total - l2_0;

l1 = l1_0 - v_rope * t;
l2 = l2_0 + v_rope * t;

%pulley positions
pulley_x = 10;
pulley_y = pulley_x * tan(theta) + 2; % Height above the slope

%x(t)
% x*cos(theta) + l1*cos(phi) = constant
x0 = 2; 
C = x0*cos(theta) + l1_0*cos(phi);
x_pos = (C - l1.*cos(phi)) ./ cos(theta);

%%graph
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; axis equal; grid on;
axis([-2 15 -2 8]);
title('Car being pulled UP Animation');

%slope drawing
slope_x_line = [-2 15];
slope_y_line = slope_x_line * tan(theta);
plot(slope_x_line, slope_y_line, 'k', 'LineWidth', 2);

%car drawing
car = plot(0, 0, 'bs', 'MarkerSize', 15, 'MarkerFaceColor', 'b');
rope_l1 = line([0 0], [0 0], 'Color', 'r', 'LineWidth', 1.5);
rope_l2 = line([0 0], [0 0], 'Color', 'r', 'LineWidth', 1.5);
pulley = plot(pulley_x, pulley_y, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [.5 .5 .5]);

%%movie magic
for k = 1:length(t)
    %car position
    cx = x_pos(k);
    cy = cx * tan(theta);
    
    set(car, 'XData', cx, 'YData', cy);
    
    set(rope_l1, 'XData', [cx, pulley_x], 'YData', [cy, pulley_y]);
    
    set(rope_l2, 'XData', [pulley_x, pulley_x], 'YData', [pulley_y, pulley_y - l2(k)*0.5]);
    
    drawnow limitrate;
    pause(1/fps);
end
