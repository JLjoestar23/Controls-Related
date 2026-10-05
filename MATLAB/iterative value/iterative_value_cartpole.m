%% Demonstrates value iteration tofind optimal swing-up control policy for a cartpole system

% discretize states
x_bins = linspace(-1, 1, 21);
dx_bins = linspace(-5, 5, 21);
theta_bins = linspace(0, 2*pi, 21);
dtheta_bins = linspace(-6*pi, 6*pi, 21);

% discretize inputs
u_bound = 10;
u = linspace(-u_bound, u_bound, 21);

% 4D discretized state-space
[x, dx, theta, dtheta] = ndgrid(x_bins, dx_bins, theta_bins, dtheta_bins);

% vectorize for computational efficiency
x_vec = x(:);
dx_vec = dx(:);
theta_vec = theta(:);
dtheta_vec = dtheta(:);

% get number of elements in each data structure
n_x = numel(x_bins);
n_dx = numel(dx_bins);
n_theta = numel(theta_bins);
n_dtheta = numel(dtheta_bins);
n_u = numel(u);
n_states = n_x * n_dx * n_theta * n_dtheta;

% construct corresponding state-space indices
[x_idx, dx_idx, theta_idx, dtheta_idx] = ndgrid(1:n_x, 1:n_dx, 1:n_theta, 1:n_dtheta);

% vector of indicies for every possible state
state_idx = sub2ind([n_x, n_dx, n_theta, n_dtheta], x_idx(:), dx_idx(:), theta_idx(:), dtheta_idx(:));
input_idx = 1:n_u;

% define goal location
goal = [0, 0, pi, 0];

%%

% define next-state matrices for every possible state-input combination
next_x_vec = zeros(n_states, n_u);
next_dx_vec = zeros(n_states, n_u);
next_theta_vec = zeros(n_states, n_u);
next_dtheta_vec = zeros(n_states, n_u);

% define additive cost matrix for every possible state-input combination
L = zeros(n_states, n_u);
%%

X = [x_vec, dx_vec, theta_vec, dtheta_vec];

next_state_idx = zeros(n_states, n_u);

%%
goal_mask = check_goal(X, goal, tol_goal);
goal_idx  = find(goal_mask);
transitions_to_goal = any(ismember(next_state_idx, goal_idx), 2);
fprintf('States that can reach goal in 1 step: %d\n', sum(transitions_to_goal));

%% calculate state transition for each state and input

tol_goal = max([x_bins(2)-x_bins(1), ...
                    dx_bins(2)-dx_bins(1), ...
                    theta_bins(2)-theta_bins(1), ...
                    dtheta_bins(2)-dtheta_bins(1)]);

Q = diag([10, 10, 1, 10]); % define quadratic state cost matrix
R = 1; % define input cost

for i=1:n_u
    % discrete change in state
    dXdt = cartpole(X, u(i));
    dt = 0.5;
    X_next = X + dXdt*dt;
    X_next(:, 3) = wrapTo2Pi(X_next(:, 3));

    x_next_idx = val_to_idx(X_next(:,1), x_bins);
    dx_next_idx = val_to_idx(X_next(:,2), dx_bins);
    theta_next_idx = val_to_idx(X_next(:,3), theta_bins);
    dtheta_next_idx = val_to_idx(X_next(:,4), dtheta_bins);

    next_state_idx(:,i) = sub2ind([n_x, n_dx, n_theta, n_dtheta], x_next_idx, dx_next_idx, theta_next_idx, dtheta_next_idx);

    % calculate cost for each state and input
    L(:, i) = cost_vec(X, u(i), goal, tol_goal);
end

%%

% initialize cost-to-go
J = zeros(n_states, 1);
tol = 0.01; % delta tolerance threshold to trigger convergence
delta = inf;
while delta > tol % iterate until convergence
    [J_next, opt_control] = min(L + 0.999*J(next_state_idx), [], 2);
    delta = max(abs(J_next - J));
    J = J_next;
    disp(delta)
end

%%
% reshape J back to 4D grid first
J_4d = reshape(J, [n_x, n_dx, n_theta, n_dtheta]);

% x-dx plane: take min over theta and dtheta dimensions (3 and 4)
J_x_dx = min(J_4d, [], 3);       % min over theta -> [n_x, n_dx, n_dtheta]
J_x_dx = min(J_x_dx, [], 3);     % min over dtheta -> [n_x, n_dx]

% theta-dtheta plane: take min over x and dx dimensions (1 and 2)
J_th_dth = min(J_4d, [], 1);     % min over x -> [n_dx, n_theta, n_dtheta]
J_th_dth = min(J_th_dth, [], 1); % min over dx -> [n_theta, n_dtheta]

% squeeze removes singleton dimensions after min operations
J_x_dx   = squeeze(min(min(J_4d, [], 3), [], 4));   % -> [n_x, n_dx]
J_th_dth = squeeze(min(min(J_4d, [], 1), [], 2));   % -> [n_theta, n_dtheta]

figure();
subplot(1,2,1);
imagesc(x_bins, dx_bins, J_x_dx');
xlabel('x (m)'); ylabel('\dot{x} (m/s)');
title('Cost-to-go: x-\dot{x} plane');
colorbar; axis xy;

subplot(1,2,2);
imagesc(theta_bins, dtheta_bins, J_th_dth');
xlabel('\theta (rad)'); ylabel('\dot{\theta} (rad/s)');
title('Cost-to-go: \theta-\dot{\theta} plane');
colorbar; axis xy;

% u_opt_x = us_vec(opt_control, 1);
% u_opt_y = us_vec(opt_control, 2);
% U_opt_x = reshape(u_opt_x, [n_y, n_x]);
% U_opt_y = reshape(u_opt_y, [n_y, n_x]);

% imagesc(xbins, ybins, J_matrix);
% hold on;
% quiver(x, y, U_opt_x, U_opt_y, 0.4, 'Color', 'k');
% axis image;
% colorbar;
% xlabel('x'); ylabel('y');
% title('Cost-to-go J* w/ Optimal Control Policy');
% hold off;

%% Functions

% define vectorized system dynamics
function dXdt = cartpole(X, u)
    mc = 0.5;
    mp = 0.2;
    l = 0.5;
    g = 9.81;
    
    x = X(:, 1);
    theta = X(:, 2);
    dx = X(:, 3);
    dtheta = X(:, 4);
    
    % vectorized 2x2 mass matrix elements
    M11 = mc + mp; % scalar
    M12 = mp*l*cos(theta); % n_states x 1
    M22 = mp*l^2; % scalar
    det_M = M11*M22 - M12.^2; % n_states x 1

    % RHS vectors
    rhs1 = u; % cart: B*u term
    rhs2 = -mp*g*l*sin(theta) + mp*l.*dtheta.^2.*sin(theta);  % pendulum: tau - C*dq

    % analytical 2x2 solve (Cramer's rule)
    ddx     = ( M22.*rhs1 - M12.*rhs2) ./ det_M;
    ddtheta = (-M12.*rhs1 + M11.*rhs2) ./ det_M;

    dXdt = [dx, dtheta, ddx, ddtheta]; % n_states x 4
end

% checks if the state is the goal
function is_goal = check_goal(X, goal, tol)
    diff    = X - goal; % n_states x 4
    is_goal = sqrt(sum(diff.^2, 2)) < tol; % n_states x 1 logical vector
end

%{
function C = cost_vec(X, u, goal, tol, Q, R)
    % quadratic state cost for all states at once
    diff = X - goal; % n_states x 4
    state_cost = sum((diff * Q) .* diff, 2); % n_states x 1  (x'Qx vectorized)

    % input cost u
    u_cost = u * R * u; % scalar

    % combine
    C = state_cost + u_cost; % n_states x 1

    % zero out goal states
    C(check_goal(X, goal, tol)) = 0;
end
%}

% simpler cost for debugging — just 1 everywhere except goal
function C = cost_vec(X, u, goal, tol_goal)
    C = ones(size(X,1), 1);
    C(check_goal(X, goal, tol_goal)) = 0;
end

% clamp to grid bounds then find nearest bin for each dimension
function idx = val_to_idx(val, bins)
    % find nearest bin index for each value in val (n_states x 1)
    val_clamped = min(max(val, bins(1)), bins(end));
    [~, idx] = min(abs(val_clamped - bins), [], 2);  % n_states x 1
end