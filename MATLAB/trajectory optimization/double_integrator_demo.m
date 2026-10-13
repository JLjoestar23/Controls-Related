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

%% Functions

function [X_sol, U_sol] = solve_quad_prog(N, A, B, x0, xf)

    nX = size(A, 1);
    nU = size(B, 2);
    nD = nX*N+nU*(N-1);
    nG = nX+nU;
    
    [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf);
    
    [H, f] = build_cost(N, nX, nU, 'quadratic');

    opts = optimoptions('quadprog', 'Display', 'none');
    z = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);
    
    initial_reshape = reshape([z; zeros(nU, 1)], nG, N);
    X_sol = initial_reshape(1:nX, :);
    U_sol = initial_reshape(nX+1:end, 1:N-1);
end

function [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf)
    % N; % number of discrete timesteps
    % nX; % state dimensions
    % nU; % input dimensions
    
    % A; % discrete-time state matrix
    % B; % discrete-time input matrix
    % x0; % system ICs

    nD = nX*N + nU*(N-1); % total number of decision variables
    nG = nX+nU; % number of decision variables per "group"
    
    % initialize equality constraint matrix and vector
    Aeq = zeros(nG*N, nD);
    beq = zeros(nG*N, 1);
    
    % setup for system initial conditions
    Aeq(1:nX, 1:nX) = eye(2);
    beq(1:nX) = x0;
    
    % dynamics equality constraints: x[k+1] = Ax[k] + Bu[k]
    for k=1:N-1
        row = nX + (k-1)*nX + (1:nX);
        x_col = (k-1)*nG + (1:nX);
        x_next_col = k*nG + (1:nX);
        u_col = x_col(end) + (1:nU);
        Aeq(row, x_col) = -A;
        Aeq(row, x_next_col) = eye(2);
        Aeq(row, u_col) = -B;
    end

    % setup for system terminal conditions
    Aeq(end-(nX-1):end, end-(nX-1):end) = eye(2);
    beq(end-(nX-1):end) = xf;
    
    % setup lower and upper bounds for states and systems
    lb = -inf(nD,1);
    ub = inf(nD,1);

    % define bounded input
    for k = 1:N-1
        u_idx = (k-1)*(nX+nU) + nX + (1:nU);
        lb(u_idx) = -1;  ub(u_idx) = 1;
    end
end

function [H, f] = build_cost(N, nX, nU, cost_type)
    if strcmpi(cost_type, 'quadratic')
        nD = nX*N+nU*(N-1);
        H = 2*eye(nD);
        f = zeros(nD, 1);
    end
end