
%% pset 6 ME 104 Triple Pendulum sim
clear; clc; close all;

%%inputs that can be edited
g_val   = 9.81;
m_val   = 1.0;
ell_val = 0.5;

%initial conditions: [a1, a2, a3, alpha1, alpha2, alpha3]
U0 = [pi/3; pi/2; 2*pi/3; 0; 0; 0];

tspan = 0:0.01:10;

syms a1 a2 a3 alpha1 alpha2 alpha3 real
syms g m ell real

%%given masses and lengths
mA = m;
mB = 2*m;
mC = 3*m;

lA = ell;
lB = 2*ell;
lC = 3*ell;

%%kinematics
% r = x e_x + y e_y
xA = lA*cos(a1);
yA = lA*sin(a1);

xB = xA + lB*cos(a2);
yB = yA + lB*sin(a2);

xC = xB + lC*cos(a3);
yC = yB + lC*sin(a3);

vA2 = lA^2*alpha1^2;

vB2 = lA^2*alpha1^2 + lB^2*alpha2^2 ...
    + 2*lA*lB*alpha1*alpha2*cos(a1-a2);

vC2 = lA^2*alpha1^2 + lB^2*alpha2^2 + lC^2*alpha3^2 ...
    + 2*lA*lB*alpha1*alpha2*cos(a1-a2) ...
    + 2*lA*lC*alpha1*alpha3*cos(a1-a3) ...
    + 2*lB*lC*alpha2*alpha3*cos(a2-a3);

%%kinetic
K = simplify(1/2*mA*vA2 + 1/2*mB*vB2 + 1/2*mC*vC2);

%%potential
U = simplify(-mA*g*xA - mB*g*xB - mC*g*xC);

%lagrangian equation
L = simplify(K - U);

%generalized coordinates
s  = [a1; a2; a3];
ds = [alpha1; alpha2; alpha3];

%solving equations
S = localEulerLagrange(s, ds, L);

%converting symbolic ang accel to handles
f_dalpha1 = matlabFunction(S.dalpha1, ...
    'Vars', {a1,a2,a3,alpha1,alpha2,alpha3,g,m,ell});
f_dalpha2 = matlabFunction(S.dalpha2, ...
    'Vars', {a1,a2,a3,alpha1,alpha2,alpha3,g,m,ell});
f_dalpha3 = matlabFunction(S.dalpha3, ...
    'Vars', {a1,a2,a3,alpha1,alpha2,alpha3,g,m,ell});

fprintf('starting ode45');


params.g = g_val;
params.m = m_val;
params.ell = ell_val;
params.f_dalpha1 = f_dalpha1;
params.f_dalpha2 = f_dalpha2;
params.f_dalpha3 = f_dalpha3;

[t,Uout] = ode45(@(t,U) triplepend_int(t,U,params), tspan, U0);

fprintf('plotting..');

%%movie magic
a1_num = Uout(:,1);
a2_num = Uout(:,2);
a3_num = Uout(:,3);

lA_num = ell_val;
lB_num = 2*ell_val;
lC_num = 3*ell_val;
Ltot = lA_num + lB_num + lC_num;

figure('Name','Triple Pend Animation');
for i = 1:length(t)
    [x, y] = get_positions(a1_num(i), a2_num(i), a3_num(i), lA_num, lB_num, lC_num);
    plot(x, y, 'o-', 'MarkerSize', 10, 'LineWidth', 1.5)
    axis equal
    axis([-Ltot Ltot -Ltot Ltot])
    grid on
    title(['Triple Pendulum, t = ', num2str(t(i), '%.2f'), ' s'])
    xlabel('y direction')
    ylabel('-x direction for plotting')
    pause(0.001)
end

%%getting five snapshots
snapshot_indices = round(linspace(1, length(t), 5));

figure('Name','Five Snapshots');
hold on
for j = 1:5
    i = snapshot_indices(j);
    [x, y] = get_positions(a1_num(i), a2_num(i), a3_num(i), lA_num, lB_num, lC_num);
    plot(x, y, 'o-', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', ['t = ', num2str(t(i), '%.2f'), ' s'])
end
axis equal
axis([-Ltot Ltot -Ltot Ltot])
grid on
legend('Location','best')
title('Five snapshots of triple pendulum motion')
xlabel('y direction')
ylabel('-x direction for plotting')
hold off

%%functions
function Udot = triplepend_int(~, U, params)
    a1 = U(1);
    a2 = U(2);
    a3 = U(3);
    alpha1 = U(4);
    alpha2 = U(5);
    alpha3 = U(6);

    da1 = alpha1;
    da2 = alpha2;
    da3 = alpha3;

    dalpha1 = params.f_dalpha1(a1,a2,a3,alpha1,alpha2,alpha3,params.g,params.m,params.ell);
    dalpha2 = params.f_dalpha2(a1,a2,a3,alpha1,alpha2,alpha3,params.g,params.m,params.ell);
    dalpha3 = params.f_dalpha3(a1,a2,a3,alpha1,alpha2,alpha3,params.g,params.m,params.ell);

    Udot = [da1; da2; da3; dalpha1; dalpha2; dalpha3];
end

function derivs = localEulerLagrange(s, ds, L)

    s = s(:);
    ds = ds(:);
    dds = str2sym("d" + string(ds));

    EQ = sym(zeros(length(s),1));

    for ii = 1:length(s)
        partial_s  = diff(L, s(ii));
        partial_ds = diff(L, ds(ii));

        partial_dt_ds = jacobian(partial_ds, [s; ds]) * [ds; dds];

        EQ(ii) = simplify(partial_dt_ds - partial_s);
    end

    derivs = solve(EQ == 0, dds);
end

function [xplot, yplot] = get_positions(a1, a2, a3, lA, lB, lC)

    xO = 0;
    yO = 0;

    xA_down = lA*cos(a1);
    yA_right = lA*sin(a1);

    xB_down = xA_down + lB*cos(a2);
    yB_right = yA_right + lB*sin(a2);

    xC_down = xB_down + lC*cos(a3);
    yC_right = yB_right + lC*sin(a3);

    xplot = [yO, yA_right, yB_right, yC_right];
    yplot = -[xO, xA_down, xB_down, xC_down];
end
