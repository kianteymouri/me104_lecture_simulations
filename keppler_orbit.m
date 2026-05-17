%%Lecture 16

function orbital_simulation_static()
    % Trajectory Equation: r(theta) = p / (1 + e * cos(theta - phi)) 

%%configuration 
    e       = 0.7;    %eccentricity
    phi_deg = 0;     %offset angle in degrees
    p       = 2.0;    %trajectory parameter


    %convert angle to radians
    phi = deg2rad(phi_deg); 

    %eval sweep
    theta = linspace(0, 2*pi, 1000);
    
    %avoid divison by zero
    denominator = 1 + e * cos(theta - phi); 
    valid_idx = denominator > 1e-4; 
    
    r = NaN(size(theta));
    r(valid_idx) = p ./ denominator(valid_idx); 

    %converting polar to cartesian
    x = r .* cos(theta);
    y = r .* sin(theta);

    %%calculating metrics
    if e == 0
        orbit_type = 'Circle'; 
        rp = p;
        ra = p;
    elseif e < 1
        orbit_type = 'Ellipse'; 
        rp = p / (1 + e); %closest
        ra = p / (1 - e); % farthest
    elseif e == 1
        orbit_type = 'Parabola'; 
        rp = p / 2; 
        ra = Inf;
    else
        orbit_type = 'Hyperbola'; 
        rp = p / (1 + e); 
        ra = Inf;
    end

    %%window readout
    fprintf('\n======================================\n');
    fprintf('ORBITAL TRAJECTORY METRICS\n');
    fprintf('======================================\n');
    fprintf('Trajectory Class : %s\n', orbit_type);
    fprintf('Parameter (p)    : %.2f\n', p); 
    fprintf('Eccentricity (e) : %.2f\n', e); 
    fprintf('Offset Angle     : %.1f°\n', phi_deg); 
    fprintf('Periapsis (r_p)  : %.2f\n', rp); 
    if isfinite(ra)
        fprintf('Apoapsis (r_a)   : %.2f\n', ra);
    else
        fprintf('Apoapsis (r_a)   : Infinity (Escape Path)\n'); 
    end
    fprintf('======================================\n');

    %%plotting
    fig = figure('Name', 'ME 104: Static Orbital Trajectory', 'Color', 'w');
    ax = axes('Parent', fig);
    grid(ax, 'on');
    axis(ax, 'equal');
    hold(ax, 'on');
    xlabel(ax, 'x (\underline{e}_x)'); 
    ylabel(ax, 'y (\underline{e}_y)'); 
    title(ax, sprintf('Orbital Shape: %s (e = %.2f)', orbit_type, e));

    %plot mass M at center
    plot(ax, 0, 0, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'Mass M'); 

    %plotting path
    plot(ax, x, y, 'r-', 'LineWidth', 2, 'DisplayName', 'Trajectory');

    % Plot the Periapsis point
    xp = rp * cos(phi);
    yp = rp * sin(phi);
    plot(ax, xp, yp, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Periapsis (r_p)');

    % Plot Apoapsis point
    if isfinite(ra)
        xa = ra * cos(phi + pi);
        ya = ra * sin(phi + pi);
        plot(ax, xa, ya, 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'DisplayName', 'Apoapsis (r_a)');
    end

    legend(ax, 'show', 'Location', 'best');

    %layout bounds
    padding = rp * 2.5;
    if isfinite(ra) && ra < 25
        padding = ra * 1.2;
    end
    xlim(ax, [-padding, padding]);
    ylim(ax, [-padding, padding]);
end
