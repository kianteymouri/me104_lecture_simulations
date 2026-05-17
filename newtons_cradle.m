%%Lecture 24

function newtons_cradle_simulation()
    clear all; close all; clc;

    %%configuration case
    Scenario = 1; % Case A: 1 ball approaching 4 stationary balls 
     %Scenario = 2; % Case B: 2 balls approaching 3 stationary balls 
    %Scenario = 3;   % Case C: 3 balls approaching 2 stationary balls
    
    num_balls = 5;  % number of balls
    R = 0.4;        %radius of balls
    L = 3.0;        %length of each string
    
    %%state arrays
    %equilibrium condition
    x_equilibrium = (0:num_balls-1) * (2 * R + 0.001);
    
    %positions and angular velocities
    theta = zeros(1, num_balls);
    omega = zeros(1, num_balls);
    
    %applying intitial conditions
    if Scenario == 1
        theta(1) = -0.4; 
    elseif Scenario == 2
        theta(1) = -0.4;  
        theta(2) = -0.4;

    else
        theta(1) = -0.4;
        theta(2) = -0.4;
        theta(3) = -0.4;  
    end

    %%plot animations
    fig = figure('Name', 'ME 104: Newton''s Cradle Simulation', 'Color', 'w', 'Position', [100, 100, 800, 400]);
    ax = axes('Parent', fig);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    
    % viewing bounds
    xlim(ax, [min(x_equilibrium) - L, max(x_equilibrium) + L]);
    ylim(ax, [-L - 0.5, 0.5]);
    title(ax, sprintf('Newton''s Cradle: Scenario %d', Scenario), 'FontSize', 12);
    
    %intiializing graph handles
    h_strings = cell(1, num_balls);
    h_balls   = cell(1, num_balls);
    
    colors = lines(num_balls); %give trace color
    for idx = 1:num_balls
        h_strings{idx} = plot(ax, [0, 0], [0, 0], 'k-', 'LineWidth', 1.5);
        h_balls{idx}   = plot(ax, 0, 0, 'o', 'MarkerSize', 24, ...
                              'MarkerFaceColor', colors(idx,:), 'MarkerEdgeColor', 'k');
    end

    %%physical parameter
    dt = 0.005;        %time step
    g = 9.81;         
    total_steps = 1200; %duration of running

    %% main loop
    for step = 1:total_steps
        
        % kinematics update
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
