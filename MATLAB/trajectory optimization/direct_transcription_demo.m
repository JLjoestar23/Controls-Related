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
N = 284;
nX = size(Ad, 1);
nU = size(Bd, 2);
x0 = [-2; 0];
xf = [0; 0];

[Aeq, beq, lb, ub] = build_constraints(N, nX, nU, Ad, Bd, x0, xf);

[H, f] = build_cost(N, nX, nU, 'quadratic');

opts = optimoptions('quadprog', 'Display', 'none');
z = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);


%% Functions

function [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf)
    % N; % number of discrete timesteps
    % nX; % state dimensions
    % nU; % input dimensions
    
    % A; % discrete-time state matrix
    % B; % discrete-time input matrix
    % x0; % system ICs

    nD = nX*N + nU*(N-1); % total number of decision variables
    
    % initialize equality constraint matrix and vector
    Aeq = zeros(nX+nX*N, nD);
    beq = zeros(nX+nX*N, 1);
    
    % setup for system initial conditions
    Aeq(1:nX, 1:nX) = eye(2);
    beq(1:nX) = x0;
    
    % dynamics equality constraints: x[k+1] = Ax[k] + Bu[k]
    for k=1:N-1
        row = nX + (k-1)*nX + (1:nX);
        x_col = (k-1)*(nX+nU) + (1:nX);
        x_next_col = k*(nX+nU) + (1:nX);
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
        H = 2*eye(N*(nX+nU));
        f = zeros(N*(nX+nU), 1);
    end
end