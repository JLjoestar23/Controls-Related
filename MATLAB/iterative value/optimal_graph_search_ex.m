%% Demonstrates value iteration to solve for a robot's optimal path to target

% define states
% 20x20 grid of possible positions
xbins = 1:20;
ybins = 1:20;
[x, y] = meshgrid(xbins, ybins);

% define inputs
[ux, uy] = meshgrid(-1:1, -1:1);

% vectorized each dimension for computational efficiency
x_vec = x(:);
y_vec = y(:);
ux_vec = ux(:);
uy_vec = uy(:);

ss_vec = [x_vec, y_vec]; % state-space vectorized
us_vec = [ux_vec, uy_vec]; % input-space vectorized

n_x = size(xbins, 2);
n_y = size(ybins, 2);
n_a = size(us_vec, 1);

% transform vector into scalar indices for reference
state_idx = sub2ind([n_y, n_x], y_vec, x_vec); 
input_idx = 1:9;

% define goal location
goal = [4, 10];

next_x_vec = zeros(n_x*n_y, n_a);
next_y_vec = zeros(n_x*n_y, n_a);
L = zeros(n_x*n_y, n_a);
for i=1:n_a
    % calculate state transition for each state and input
    next_x_vec(:, i) = min(max(x_vec + us_vec(i, 1), 1), 20);
    next_y_vec(:, i) = min(max(y_vec + us_vec(i, 2), 1), 20);
    
    % calculate cost for each state and input
    L(:, i) = cost_vec(ss_vec, us_vec(i, :), goal);
end
next_state_idx = sub2ind([n_y, n_x], next_y_vec, next_x_vec);

% initialize cost-to-go
J = zeros(n_x*n_y, 1);
tol = 0.01; % delta tolerance threshold to trigger convergence
delta = inf;
while delta > tol % iterate until convergence
    [J_next, opt_control] = min(L + J(next_state_idx), [], 2);
    delta = max(abs(J_next - J));
    J = J_next;
end
J_matrix = reshape(J, [n_y, n_x]);
u_opt_x = us_vec(opt_control, 1);
u_opt_y = us_vec(opt_control, 2);
U_opt_x = reshape(u_opt_x, [n_y, n_x]);
U_opt_y = reshape(u_opt_y, [n_y, n_x]);

imagesc(xbins, ybins, J_matrix);
hold on;
quiver(x, y, U_opt_x, U_opt_y, 0.4, 'Color', 'k');
axis image;
colorbar;
xlabel('x'); ylabel('y');
title('Cost-to-go J* w/ Optimal Control Policy');
hold off;

%% Functions

% define obstacles
function is_obstacle = check_obstacle(X)
    is_obstacle = zeros(size(X, 1), 1);
    for i=1:size(X, 1)
        if X(i, 1) > 6 && X(i, 1) < 10 && X(i, 2) > 4 && X(i, 2) < 16
            is_obstacle(i) = true;
        end
    end
end

% define cost based on state and action
function C = cost_vec(X, u, goal)
    C = zeros(size(X, 1), 1);
    for i=1:size(X, 1)
        if isequal(X(i,:), goal) % goal has 0 cost
            C(i) = 0;
        elseif check_obstacle(X(i,:)) % obstacles have high cost
            C(i) = 10;
        else % normal states have 1 cost
            C(i) = 1;
        end
        % heavy cost for diagonal movements
        % strictly vertical and horizontal moves
        u_cost = norm(u,1);
        if u_cost > 1
            u_cost = 10;
        end
        C(i) = C(i) + u_cost; % combine state and input cost for total
    end
end