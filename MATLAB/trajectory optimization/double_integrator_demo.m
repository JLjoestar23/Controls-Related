%% Direct Transcription Demo for a Double Integrator System
% Demonstrates the implementation and solving for the optimal trajectory
% and control policy for an arbitrary double integrator system

% Continuous time dynamics
A = [0 1; 0 0];
B = [0; 1];
C = eye(2);
D = 0;
sys_ct = ss(A, B, C, D);

% Discrete time dynamics
% Can either discretize using MATLAB's c2d() or a first-order Euler 
% approximation
Ts = 0.01; % sampling period
discretize_method = 'euler'; % change dicretization method
if strcmpi(discretize_method, 'c2d')
    sys_dt = c2d(sys_ct, Ts);
    Ad = sys_dt.A;
    Bd = sys_dt.B;
elseif strcmpi(discretize_method, 'euler')
    Ad = eye(2) + A*Ts;
    Bd = B*Ts;
end

% Establish system initial and terminal conditions
x0 = [-2; 0];
xf = [10; 0];

% Find minimum-time N via bisection search
% Ts is fixed, since if Ts was a decision variable for a strictly
% minimum-time problem, system dynamics would be nonlinear, resulting in a 
% non-convex problem, to which quadprog cannot be used. Instead, a search 
% for the minimum amount of timesteps required for a feasible solution
% produces a minimum-time optimal trajectory and control policy for a
% quadratic cost function.
N_lo = 2;
N_hi = 1000;
while N_hi - N_lo > 1
    N_mid = floor((N_lo + N_hi)/2);
    feasible = is_feasible(N_mid, Ad, Bd, x0, xf);
    if feasible
        N_hi = N_mid;
    else
        N_lo = N_mid;
    end
end
N_min = N_hi; % minimum feasible N
fprintf('N_min = %d, t_f = %.4f s\n', N_min, (N_min-1)*Ts);

% Re-solve at N_min as a quadratic program
N = N_min;
[X_sol, U_sol, fval, exitflag, output] = solve_quadprog(N, Ad, Bd, x0, xf);

% Plotting results
close all;

% Create time arrays based on Ts
t_X = (0:N-1)*Ts;
t_U = (0:N-2)*Ts;

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