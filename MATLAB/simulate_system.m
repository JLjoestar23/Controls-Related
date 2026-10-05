%% Double Pendulum
% Parameters
system_params.m1 = 2; 
system_params.m2 = 2;
system_params.l1 = 1;
system_params.l2 = 1;
system_params.g = 9.81;

% Simulate
tspan = [0 20];
X0    = [pi/2; pi/2; 0; 0];
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, X] = ode45(@(t,X) double_pendulum(t, X, system_params), tspan, X0, opts);

%% Animate Double Pendulum
double_pendulum_anim(X, t, tspan, system_params);

%% Single Pendulum
% Parameters
system_params.m = 2; 
system_params.l = 1;
system_params.g = 9.81;

% Simulate
tspan = [0 20];
X0    = [0.1; 0];
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, X] = ode45(@(t,X) single_pendulum(t, X, system_params), tspan, X0, opts);

%% Animate Double Pendulum
single_pendulum_anim(X, t, tspan, system_params);

%% Cartpole Uncontrolled
% Parameters
system_params.mc = 2;
system_params.mp = 1;
system_params.l = 1;
system_params.g = 9.81;

% Simulate
tspan = [0 5];
X0    = [0; pi/8; 0; 0]; % Initial state for cartpole
opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, X] = ode45(@(t,X) cartpole(t, X, system_params), tspan, X0, opts);

%% Animate Cartpole
cartpole_anim(X, t, tspan, system_params);

%% Cartpole LQR Design
mc = 0.5;
mp = 0.2;
l = 0.5;
g = 9.81;

system_params.mc = mc;
system_params.mp = mp;
system_params.l = l;
system_params.g = g;

A = [0      0                   1     0;
     0      0                   0     1;
     0      mp*g/mc             0     0;
     0      (mc+mp)*g/(mc*l)    0     0];
B = [0; 0; 1/mc; 1/(mc*l)];
C = eye(4);
D = 0;

cartpole_sys = ss(A, B, C, D);

x_max = 0.2; % m
dx_max = 0.5; % m/s
theta_max = 0.1; % rad
dtheta_max = 1; % rad/s

Q = diag([1/x_max^2; 1/theta_max^2; 1/dx_max^2; 1/dtheta_max^2]);
R = 1/(40^2);

K = lqr(cartpole_sys, Q, R);
K = round(K, 4);
% system_params.K = zeros(1, 4);
system_params.K = K;

%% Simulate
close all;
clc;

tspan = [0 10];
X0    = [1; -0.1; 0; 0]; % Initial state for cartpole
% opts  = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, X] = ode45(@(t,X) cartpole(t, X, system_params), tspan, X0);

figure;
plot(t, X);
hold on;
legend('x', 'theta', 'dx', 'dtheta');
grid on;
hold off;

dtheta = X(:, 4);
theta = X(:, 2);

E = 0.5*mp*l^2*dtheta.^2 - mp*g*l*cos(theta);
E_des = mp*g*l;
E_err = E - E_des;

figure;
plot(t, E);
hold on;
plot(t, E_des*ones(size(t)));
plot(t, E_err);
hold off;

cartpole_anim(X, t, tspan, system_params);
