%%Lecture 28

function rolling_wheels_rigid_bar_universal()
    clear all; close all; clc;

    %%configuration parameter
    R1 = 1.0;        % wheel radius 1
    R2 = 1.5;        %wheel radius 2
    D0 = 4.0;        % initial horizontal distance
    
    %kinematics for wheel 1
    omega1_init = 1.5;  % Initial angular velocity of wheel 1 (rad/s)
    alpha1      = 0.2;  % cnst angular acceleration of wheel 1 (rad/s^2)

    %% time parameters
    dt = 0.01;          
    t_max = 5.0;        
    time = 0:dt:t_max;
    
    %%intialize geometrey
    L_AB = sqrt(D0^2 + ( (R2 + R2/2) - (R1 + R1/2) )^2);
    
    % track positions over time
    xc1 = 0;           
    theta1 = 0;         % Angle of wheel 1
    omega1 = omega1_init;

    %% setting up graphics
    fig = figure('Name', 'ME 104: Perfectly Rigid Bar Simulation', 'Color', 'w', 'Position', [100, 100, 900, 450]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    
    xlabel(ax, 'x (\underline{e}_x)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ax, 'y (\underline{e}_y)', 'Interpreter', 'latex', 'FontSize', 12);

    %%movie magic
    for t = time
        %%updating wheel 1
        vc1 = -omega1 * R1;           %cneter 1 velo
        xc1 = xc1 + vc1 * dt;         % moving center 1
        theta1 = theta1 + omega1 * dt; %rotating wheel 1
        
        C1 = [xc1, R1];
        % position of Pin A
        A = C1 + [-(R1/2)*sin(theta1), (R1/2)*cos(theta1)];

        %% rigid link kinematics
        theta2_guess = theta1 * (R1/R2); 
        
        options = optimset('Display', 'off', 'TolX', 1e-6);
        theta2 = fzero(@(t2) distance_error_root(t2, A, xc1, D0, R2, L_AB), theta2_guess, options);
        
        %reconstructing wheel 2 center
        xc2 = D0 - theta2 * R2; 
        C2 = [xc2, R2];
        B = C2 + [-(R2/2)*sin(theta2), (R2/2)*cos(theta2)];

        %%metrics calc
        current_bar_length = norm(B - A);

        %%drawing graphics
        cla(ax); 
        plot(ax, [-20, 25], [0, 0], 'k-', 'LineWidth', 2); % Floor
        
        % drawing wheels
        draw_wheel(ax, C1, R1, theta1, [0.2 0.6 0.8]);
        draw_wheel(ax, C2, R2, theta2, [0.8 0.4 0.2]);
        
        %draw bar
        plot(ax, [A(1), B(1)], [A(2), B(2)], 'g-', 'LineWidth', 4);
        
        % draw pins
        plot(ax, A(1), A(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
        text(ax, A(1)-0.2, A(2)+0.3, 'A', 'FontWeight', 'bold');
        plot(ax, B(1), B(2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
        text(ax, B(1)-0.2, B(2)+0.3, 'B', 'FontWeight', 'bold');

        %centering camera
        xlim(ax, [xc1 - 3*R1, xc2 + 3*R2]);
        ylim(ax, [-0.5, max(R1, R2) * 2.8]);
        
        title_str = sprintf('Time: %.2fs | Bar Length: %.4f (Constant!)', t, current_bar_length);
        title(ax, title_str, 'FontSize', 11);
        
        drawnow;
        
        %advance driver state
        omega1 = omega1 + alpha1 * dt; 
    end
end


function error_val = distance_error_root(theta2, A, xc1, D0, R2, L_AB)
    xc2 = D0 - theta2 * R2;
    
    B = [xc2, R2] + [-(R2/2)*sin(theta2), (R2/2)*cos(theta2)];
    
    error_val = norm(B - A) - L_AB;
end

%%rendering elements
function draw_wheel(ax, center, radius, angle, color)
    theta_circle = linspace(0, 2*pi, 100);
    x_rim = center(1) + radius * cos(theta_circle);
    y_rim = center(2) + radius * sin(theta_circle);
    plot(ax, x_rim, y_rim, 'Color', color, 'LineWidth', 2.5);
    
    num_spokes = 4;
    for i = 1:num_spokes
        spoke_angle = angle + (i - 1) * (pi / 2) + pi/2;
        x_spoke = center(1) + radius * cos(spoke_angle);
        y_spoke = center(2) + radius * sin(spoke_angle);
        plot(ax, [center(1), x_spoke], [center(2), y_spoke], 'k-', 'LineWidth', 1);
    end
    plot(ax, center(1), center(2), 'k.', 'MarkerSize', 10);
end
