%% Demonstrates the use of value iteration to find the optimal swing-up 
% control policy for a torque-controlled pendulum

% discretize states
theta_bins = linspace(0, 2*pi, 301);
dtheta_bins = linspace(-6*pi, 6*pi, 301);

% discretize input-space
u_bound = 10;
u = linspace(-u_bound, u_bound, 301);

% discretized state-space
[theta, dtheta] = ndgrid(theta_bins, dtheta_bins);

% vectorize for computational steps later on
theta_vec = theta(:);
dtheta_vec = dtheta(:);

% get number of elements in each data structure
n_theta = numel(theta_bins);
n_dtheta = numel(dtheta_bins);
n_u = numel(u);
n_states = n_theta * n_dtheta;

% construct corresponding state-space indices
[theta_idx, dtheta_idx] = ndgrid(1:n_theta, 1:n_dtheta);

% vector of indices for every possible state and input
% state_idx = sub2ind([n_theta, n_dtheta], theta_idx(:), dtheta_idx(:));
state_idx = sub2ind([n_theta, n_dtheta], theta_idx(:), dtheta_idx(:));
input_idx = 1:n_u;

% define goal state
% in this example, the goal state is the pendulum's upright, unstable
% equillibirum
goal = [pi, 0];

%% calculate state transition for each state and input

% these next-state matrices will store every possible state at the next
% timestep for every possible input
next_theta_vec = zeros(n_states, n_u);
next_dtheta_vec = zeros(n_states, n_u);

% define the additive cost matrix for every possible state-input 
% combination
L = zeros(n_states, n_u);

% matrix containing every possible current state
X = [theta_vec, dtheta_vec];

% this vector is to store the corresponding indices of the next-states
next_state_idx = zeros(n_states, n_u);

% set tolerance to max difference between bins
tol_goal = max([theta_bins(2)-theta_bins(1), dtheta_bins(2)-dtheta_bins(1)]);

Q = diag([1, 0]); % define quadratic state cost matrix
R = 0.005; % define input cost

for i=1:n_u
    % discrete change in state
    dXdt = t_pend(X, u(i));
    dt = 0.1; % 0.1s timestep
    X_next = X + dXdt*dt;

    % angle wrap to 2pi
    X_next(:, 1) = wrapTo2Pi(X_next(:, 1));
    
    % identify out-of-bounds before clamping dtheta
    out_of_bounds = X_next(:,2) < dtheta_bins(1) | ...
                    X_next(:,2) > dtheta_bins(end);

    % clamp to keep indices valid
    X_next(:,2) = min(max(X_next(:,2), dtheta_bins(1)), dtheta_bins(end));
    
    theta_next_idx = val_to_idx(X_next(:,1), theta_bins);
    dtheta_next_idx = val_to_idx(X_next(:,2), dtheta_bins);
    
    next_state_idx(:,i) = sub2ind([n_theta, n_dtheta], theta_next_idx, dtheta_next_idx);
    
    % calculate cost for each state and input
    % L(:, i) = cost_vec(X, u(i), goal, tol_goal);
    L(:, i) = quadratic_cost_vec(X, u(i), goal, tol_goal, Q, R);
    % L(out_of_bounds, i) = L(out_of_bounds, i) + 100;
end

%%

% initialize cost-to-go matrix
J = zeros(n_states, 1);
tol = 0.05; % delta tolerance threshold to trigger convergence
delta = inf;
iters = 0;
n_chars = 0;
y = 1; % discount

% iterative value loop to solve for J*
while delta > tol % iterate until convergence
    [J_next, opt_control] = min(L + y*J(next_state_idx), [], 2);
    delta = max(abs(J_next - J));
    J = J_next;
    iters = iters + 1;

    % erase previous line
    fprintf(repmat('\b', 1, n_chars));

    % print new status and record length
    msg = sprintf('Iteration: %d  |  Delta: %.3f', iters, delta);
    n_chars = length(msg);
    fprintf(msg);
end
fprintf('\nValue iteration converged in %d iterations\n', iters);

%% Plotting results

% plot optimal time-to-go based on state
J_matrix = reshape(J, [n_theta, n_dtheta]);

figure;
subplot(1, 2, 1);
imagesc(theta_bins, dtheta_bins, J_matrix');
hold on;
axis xy;
axis square;
colorbar;
xlabel('$$\theta$$', 'Interpreter', 'latex');
ylabel('$$\dot{\theta}$$', 'Interpreter', 'latex');
title('Cost-To-Go $J^*$', 'Interpreter', 'latex');
hold off;

% plot optimal control policy based on state
opt_control_idx_matrix = reshape(opt_control, [n_theta, n_dtheta]);
u_opt_matrix = u(opt_control_idx_matrix);

subplot(1, 2, 2);
imagesc(theta_bins, dtheta_bins, u_opt_matrix');
hold on;
axis xy;
axis square;
colorbar;
xlabel('$$\theta$$', 'Interpreter', 'latex');
ylabel('$$\dot{\theta}$$', 'Interpreter', 'latex');
title('Optimal Control Policy $\pi^*$', 'Interpreter', 'latex');
hold off;

%% Functions

% define vectorized system dynamics
function dXdt = t_pend(X, u)
    m = 2;
    l = 1;
    g = 9.81;

    theta  = X(:, 1);
    dtheta  = X(:, 2);

    M = m*l^2;
    
    F_g = m*g*l*sin(theta);
    
    b = 0.1;
    
    ddtheta = (-(F_g + b*dtheta) + u) / M;

    dXdt   = [dtheta, ddtheta]; % nx2, each column is a vector of states
end

% checks if the state is the goal
function is_goal = check_goal(X, goal, tol)
    % wrap angular error to [-pi, pi]
    ang_err  = wrapToPi(X(:,1) - goal(1));
    vel_err  = X(:,2) - goal(2);
    is_goal  = sqrt(ang_err.^2 + vel_err.^2) < tol;
end

function C = quadratic_cost_vec(X, u, goal, tol_goal, Q, R)
    % shift theta so goal is at origin
    diff       = [wrapToPi(X(:,1) - goal(1)), X(:,2) - goal(2)];
    % state_cost = 2 * sum(diff.^2, 2); % 2*x'x
    state_cost = sum((diff * Q) .* diff, 2); % x'Qx
    u_cost     = R*u^2; % u'Ru
    C          = state_cost + u_cost;
    C(check_goal(X, goal, tol_goal)) = 0;
end

% simpler cost for debugging — just 1 everywhere except goal
function C = cost_vec(X, u, goal, tol_goal)
    C = 5*ones(size(X,1), 1);
    C(check_goal(X, goal, tol_goal)) = 0;
end

% clamp to grid bounds then find nearest bin for each dimension
function idx = val_to_idx(val, bins)
    val_clamped = min(max(val, bins(1)), bins(end));
    idx = round(interp1(bins, 1:numel(bins), val_clamped));
end