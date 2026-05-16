function newtons_cradle_simulation()
    % ME 104: Newton's Cradle Discrete Collision Simulation
    % Based on Lec 24 binary impact sequence (e = 1, identical masses) [cite: 105, 112]
    clear all; close all; clc;

    %% --- Configuration Zone ---
    % Choose your lecture experiment scenario:
    Scenario = 1; % Case A: 1 ball approaching 4 stationary balls [cite: 106]
     %Scenario = 2; % Case B: 2 balls approaching 3 stationary balls [cite: 128]
    %Scenario = 3;   % Case C: 3 balls approaching 2 stationary balls
    
    num_balls = 5;  % Total number of identical balls [cite: 105]
    R = 0.4;        % Radius of each ball
    L = 3.0;        % Length of each supporting pendulum string
    
    %% --- Initialize State Arrays ---
    % Equilibrium X coordinates where balls rest side-by-side with tiny gaps [cite: 108]
    x_equilibrium = (0:num_balls-1) * (2 * R + 0.001);
    
    % Positions (angles in radians) and angular velocities (rad/s)
    theta = zeros(1, num_balls);
    omega = zeros(1, num_balls);
    
    % Apply initial conditions based on the selected lecture scenario
    if Scenario == 1
        % 1 ball approaching rightward [cite: 106]
        theta(1) = -0.4;  % Leftmost ball pulled back
    elseif Scenario == 2
        % 2 balls approaching rightward [cite: 128]
        theta(1) = -0.4;  
        theta(2) = -0.4;  % First two balls pulled back together

    else
        % Scenario 3: 3 balls approaching rightward
        theta(1) = -0.4;
        theta(2) = -0.4;
        theta(3) = -0.4;  % First three balls pulled back together
    end

    %% --- Setup Plot Animation Window ---
    fig = figure('Name', 'ME 104: Newton''s Cradle Simulation', 'Color', 'w', 'Position', [100, 100, 800, 400]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    
    % Set stable viewing bounds around the cradle geometry
    xlim(ax, [min(x_equilibrium) - L, max(x_equilibrium) + L]);
    ylim(ax, [-L - 0.5, 0.5]);
    title(ax, sprintf('Newton''s Cradle: Scenario %d', Scenario), 'FontSize', 12);
    
    % Initialize graphical handles for strings and ball masses
    h_strings = cell(1, num_balls);
    h_balls   = cell(1, num_balls);
    
    colors = lines(num_balls); % Give each ball a distinct trace color
    for idx = 1:num_balls
        h_strings{idx} = plot(ax, [0, 0], [0, 0], 'k-', 'LineWidth', 1.5);
        h_balls{idx}   = plot(ax, 0, 0, 'o', 'MarkerSize', 24, ...
                              'MarkerFaceColor', colors(idx,:), 'MarkerEdgeColor', 'k');
    end

    %% --- Physics Engine Parameters ---
    dt = 0.005;         % Simulation time step
    g = 9.81;           % Gravitational acceleration constant
    total_steps = 1200; % Duration of active running window

    %% --- Main Simulation Loop ---
    for step = 1:total_steps
        
        % 1. Kinematics Update (Standard Nonlinear Pendulum integration)
        for i = 1:num_balls
            alpha = -(g / L) * sin(theta(i)); % Angular acceleration
            omega(i) = omega(i) + alpha * dt;   % Update angular velocity
            theta(i) = theta(i) + omega(i) * dt; % Update angular position
        end
        
        % 2. Discrete Binary Collision Resolution Loop [cite: 132]
        % We check multiple passes per frame to capture chained sequence reactions correctly
        collision_detected = true;
        while collision_detected
            collision_detected = false;
            
            for i = 1:(num_balls - 1)
                % Compute current spatial X positions to verify contact boundary
                x_left  = x_equilibrium(i)   + L * sin(theta(i));
                x_right = x_equilibrium(i+1) + L * sin(theta(i+1));
                
                % Check if adjacent balls overlap or strike each other
                if x_left + R >= x_right - R
                    
                    % Is the left ball physically moving toward the right ball?
                    v_left  = omega(i)   * L;
                    v_right = omega(i+1) * L;
                    
                    if v_left > v_right
                        % Perfect elastic momentum transfer (e = 1) for identical masses: [cite: 107, 112]
                        % Velocities are completely swapped [cite: 113, 114]
                        v_left_new  = v_right;
                        v_right_new = v_left;
                        
                        % Translate linear velocities back into angular domain
                        omega(i)   = v_left_new / L;
                        omega(i+1) = v_right_new / L;
                        
                        % Prevent overlapping by separating them exactly to contact distance
                        midpoint = (x_left + x_right) / 2;
                        theta(i)   = asin(((midpoint - R) - x_equilibrium(i)) / L);
                        theta(i+1) = asin(((midpoint + R) - x_equilibrium(i+1)) / L);
                        
                        collision_detected = true; % Re-verify cascade interactions [cite: 132]
                    end
                end
            end
        end
        
        % 3. Graphics Rendering Update
        for i = 1:num_balls
            % Compute instantaneous positions relative to each anchor point
            x_pos = x_equilibrium(i) + L * sin(theta(i));
            y_pos = -L * cos(theta(i));
            
            % Update line coordinate endpoints (pendulum strings)
            set(h_strings{i}, 'XData', [x_equilibrium(i), x_pos], 'YData', [0, y_pos]);
            
            % Update marker positions (spherical ball objects)
            set(h_balls{i}, 'XData', x_pos, 'YData', y_pos);
        end
        
        drawnow;
        pause(0.001); % Slow execution down slightly for smooth viewing
    end
end
