%% pset 6 ME 104 rolling cylinders in half pipe
clear; clc; close all;

%% editable inputs
g_val = 9.81;
m_val = 1.0;
R_val = 2.0;
r_val = 0.2;
k_val = 80.0;
ell_val = 1.0;

% initial conditions: [alpha, beta, alphadot, betadot]
U0 = [0.45; 1.35; 0; 0];

tspan = 0:0.01:10;

%% symbolic variables
syms alpha beta alphadot betadot real
syms alphaddot betaddot real
syms g m R r k ell real

%% generalized coordinates
q  = [alpha; beta];
dq = [alphadot; betadot];
ddq = [alphaddot; betaddot];

%% positions
xA = R*cos(alpha);
yA = R*sin(alpha);

xB = R*cos(beta);
yB = R*sin(beta);

%% velocities squared
vA2 = R^2*alphadot^2;
vB2 = R^2*betadot^2;

%% rolling without slipping
omega1 = R/r * alphadot;
omega2 = R/r * betadot;

%% solid cylinder inertia
I = 1/2*m*r^2;

%% kinetic energy
K = simplify( ...
    1/2*m*vA2 + 1/2*I*omega1^2 + ...
    1/2*m*vB2 + 1/2*I*omega2^2 );

%% spring length
spring_len = simplify(sqrt((xB - xA)^2 + (yB - yA)^2));

%% potential energy
Ug = simplify(-m*g*xA - m*g*xB);
Us = simplify(1/2*k*(spring_len - ell)^2);

U = simplify(Ug + Us);

%% Lagrangian
L = simplify(K - U);

%% Euler-Lagrange equations
S = localEulerLagrange(q, dq, ddq, L);

f_alphaddot = matlabFunction(S.alphaddot, ...
    'Vars', {alpha,beta,alphadot,betadot,g,m,R,r,k,ell});

f_betaddot = matlabFunction(S.betaddot, ...
    'Vars', {alpha,beta,alphadot,betadot,g,m,R,r,k,ell});

%% parameters
params.g = g_val;
params.m = m_val;
params.R = R_val;
params.r = r_val;
params.k = k_val;
params.ell = ell_val;
params.f_alphaddot = f_alphaddot;
params.f_betaddot = f_betaddot;

%% solve ODE
[t,Uout] = ode45(@(t,U) halfpipe_int(t,U,params), tspan, U0);

%% extract motion
alpha_num = Uout(:,1);
beta_num  = Uout(:,2);

%% animation
figure('Name','Rolling Cylinders in Half Pipe');

for i = 1:length(t)

    [xA_plot, yA_plot, xB_plot, yB_plot, xpipe, ypipe] = ...
        get_positions(alpha_num(i), beta_num(i), R_val);

    plot(ypipe, -xpipe, 'k', 'LineWidth', 2)
    hold on

    plot(yA_plot, -xA_plot, 'ro', 'MarkerSize', 18, 'LineWidth', 2)
    plot(yB_plot, -xB_plot, 'bo', 'MarkerSize', 18, 'LineWidth', 2)

    plot([yA_plot yB_plot], -[xA_plot xB_plot], 'g-', 'LineWidth', 2)

    hold off
    axis equal
    axis([-R_val R_val -R_val 0.5])
    grid on

    title(['Rolling cylinders, t = ', num2str(t(i),'%.2f'), ' s'])
    xlabel('e_y direction')
    ylabel('-e_x direction')

    pause(0.001)
end

%% snapshots
snapshot_indices = round(linspace(1, length(t), 6));

figure('Name','Snapshots');
hold on

[xA_plot, yA_plot, xB_plot, yB_plot, xpipe, ypipe] = ...
    get_positions(alpha_num(1), beta_num(1), R_val);

% plot half pipe
h_pipe = plot(ypipe, -xpipe, 'k', 'LineWidth', 2, ...
    'DisplayName', 'Half-pipe');

for j = 1:length(snapshot_indices)
    i = snapshot_indices(j);

    [xA_plot, yA_plot, xB_plot, yB_plot, ~, ~] = ...
        get_positions(alpha_num(i), beta_num(i), R_val);

    time_label = ['t = ', num2str(t(i), '%.2f'), ' s'];

    hA = plot(yA_plot, -xA_plot, 'ro', ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'DisplayName', ['Cylinder A, ', time_label]);

    hB = plot(yB_plot, -xB_plot, 'bo', ...
        'MarkerSize', 12, 'LineWidth', 1.5, ...
        'DisplayName', ['Cylinder B, ', time_label]);

    hS = plot([yA_plot yB_plot], -[xA_plot xB_plot], ...
        'LineWidth', 1.5, ...
        'DisplayName', ['Spring, ', time_label]);
end

axis equal
axis([-R_val R_val -R_val 0.5])
grid on
title('Snapshots of rolling cylinder motion')
xlabel('e_y direction')
ylabel('-e_x direction')

legend('Location','bestoutside')

hold off

%% functions

function Udot = halfpipe_int(~, U, params)

    alpha = U(1);
    beta = U(2);
    alphadot = U(3);
    betadot = U(4);

    dalpha = alphadot;
    dbeta = betadot;

    dalphadot = params.f_alphaddot(alpha,beta,alphadot,betadot, ...
        params.g,params.m,params.R,params.r,params.k,params.ell);

    dbetadot = params.f_betaddot(alpha,beta,alphadot,betadot, ...
        params.g,params.m,params.R,params.r,params.k,params.ell);

    Udot = [dalpha; dbeta; dalphadot; dbetadot];
end

function derivs = localEulerLagrange(q, dq, ddq, L)

    q = q(:);
    dq = dq(:);
    ddq = ddq(:);

    EQ = sym(zeros(length(q),1));

    for ii = 1:length(q)

        partial_q  = simplify(diff(L, q(ii)));
        partial_dq = simplify(diff(L, dq(ii)));

        partial_dt_dq = simplify(jacobian(partial_dq, [q; dq]) * [dq; ddq]);

        EQ(ii) = simplify(partial_dt_dq - partial_q);

        fprintf('\nEuler-Lagrange equation %d:\n', ii);
        disp(latex(EQ(ii)))
    end

    derivs = solve(EQ == 0, ddq);

    fprintf('\nSolution for alphaddot:\n');
    disp(latex(simplify(derivs.alphaddot)))

    fprintf('\nSolution for betaddot:\n');
    disp(latex(simplify(derivs.betaddot)))
end

function [xA, yA, xB, yB, xpipe, ypipe] = get_positions(alpha, beta, R)

    xA = R*cos(alpha);
    yA = R*sin(alpha);

    xB = R*cos(beta);
    yB = R*sin(beta);

    theta = linspace(-pi/2, pi/2, 300);

    xpipe = R*cos(theta);
    ypipe = R*sin(theta);
end
