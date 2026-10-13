%% Direct Transcription Demo for a Double Integrator System
% Demonstrates the implementation and solving for the optimal trajectory
% and control policy for a arbitrary double integrator system

%% Establish System Dynamics

% continuous time dynamics
A = [0 1; 0 0];
B = [0; 1];
C = eye(2);
D = 0;
sys_ct = ss(A, B, C, D);

% discrete time dynamics
% using matlabs c2d
% Ts = 0.01; % s
% sys_dt = c2d(sys_ct, Ts);
% using first order Euler approximation
Ts = 0.01;
Ad = eye(2) + Ts * [0, 1; 0, 0];
Bd = Ts * [0; 1];

% Problem dimensions and boundary conditions
N = 300;
% nX = size(Ad, 1);
% nU = size(Bd, 2);
% nD = nX*N+nU*(N-1);
% nG = nX+nU;

x0 = [-2; 1];
xf = [0; 0];

[X_sol, U_sol] = solve_quad_prog(N, Ad, Bd, x0, xf);

t_X = (0:N-1)*Ts;
t_U = (0:N-2)*Ts;

figure;
plot(X_sol(1,:), X_sol(2,:), 'LineWidth', 2);
xlabel('x'); ylabel('xdot');
axis square; axis equal;
grid on; grid minor;
title('Phase Portrait');

figure;
subplot(2,1,1);
plot(t_X, X_sol(1,:), t_X, X_sol(2,:), 'LineWidth', 2);
xlabel('t (s)'); ylabel('state');
grid on;
legend('x','xdot');
title('Optimal State and Input Trajectories');

subplot(2,1,2);
plot(t_U, U_sol, 'LineWidth', 2);
xlabel('t (s)'); ylabel('u');
grid on;