%% Demonstrates the use of value iteration to find the optimal swing-up 
% control policy for a torque-controlled pendulum
% -- GPU-accelerated version (Parallel Computing Toolbox) --

%% Check GPU availability
useGPU = gpuDeviceCount("available") > 0;
if useGPU
    gpuDevice(1); % select first available GPU, prints device info
    fprintf('Using GPU: %s\n', gpuDevice().Name);
    device_str_precompute = 'GPU';
else
    warning('No GPU detected / Parallel Computing Toolbox unavailable — falling back to CPU.');
    device_str_precompute = 'CPU';
end

% discretize states
theta_bins = linspace(0, 2*pi, 601);
dtheta_bins = linspace(-4*pi, 4*pi, 601);

% discretize input-space
u_bound = 6;
u = linspace(-u_bound, u_bound, 601);

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
state_idx = sub2ind([n_theta, n_dtheta], theta_idx(:), dtheta_idx(:));
input_idx = 1:n_u;

% define goal state
% in this example, the goal state is the pendulum's upright, unstable
% equillibirum
goal = [pi, 0];

%% calculate state transition for each state and input
% At 601^3 the precompute arrays are now large enough (n_states x n_u =
% 361,201 x 601 ≈ 217M elements per matrix) that this loop is worth
% moving to the GPU too. Everything inside the loop body — t_pend,
% wrapTo2Pi_manual/wrapToPi_manual, val_to_idx_uniform, sub2ind,
% quadratic_cost_vec — is elementwise/trig or index arithmetic, all of
% which is gpuArray-safe, so no function bodies need to change.
%
% Memory note: L and next_state_idx are each ~217M elements. As double
% (8 bytes) that's ~1.7 GB per array on the GPU; as single (4 bytes,
% precisionClass below) it's ~0.87 GB. Check gpuDevice().AvailableMemory
% before running, and drop to 'single' if you're tight on VRAM.

precisionClass = 'double'; % switch to 'single' if you hit out-of-memory

if useGPU
    castfn = @(x) gpuArray(cast(x, precisionClass));
else
    castfn = @(x) cast(x, precisionClass);
end

% matrix containing every possible current state — pushed to GPU once;
% every op derived from X (t_pend, wraps, cost, indexing) then runs
% on-device automatically because X is a gpuArray.
X = castfn([theta_vec, dtheta_vec]);

% preallocate directly on the GPU (as gpuArray from the start, not
% zeros() then transferred) so no host-to-device copy happens per-input
L = castfn(zeros(n_states, n_u));
next_state_idx = castfn(zeros(n_states, n_u));

% set tolerance to max difference between bins
tol_goal = max([theta_bins(2)-theta_bins(1), dtheta_bins(2)-dtheta_bins(1)]);

Q = diag([1, 0]); % define quadratic state cost matrix
R = 0.1; % define input cost

tic;
for i=1:n_u
    % discrete change in state — u(i) is a plain host scalar; MATLAB
    % implicitly expands it against the gpuArray X with no extra transfer
    dXdt = t_pend(X, u(i));
    dt = 0.1; % 0.1s timestep
    X_next = X + dXdt*dt;

    % angle wrap to 2pi (manual — avoids Mapping Toolbox dependency,
    % and runs on-device since X_next is a gpuArray)
    X_next(:, 1) = wrapTo2Pi_manual(X_next(:, 1));
    
    % identify out-of-bounds before clamping dtheta
    out_of_bounds = X_next(:,2) < dtheta_bins(1) | ...
                    X_next(:,2) > dtheta_bins(end);

    % clamp to keep indices valid
    X_next(:,2) = min(max(X_next(:,2), dtheta_bins(1)), dtheta_bins(end));
    
    theta_next_idx = val_to_idx_uniform(X_next(:,1), theta_bins);
    dtheta_next_idx = val_to_idx_uniform(X_next(:,2), dtheta_bins);
    
    next_state_idx(:,i) = sub2ind([n_theta, n_dtheta], theta_next_idx, dtheta_next_idx);
    
    % calculate cost for each state and input
    L(:, i) = quadratic_cost_vec(X, u(i), goal, tol_goal, Q, R);
end
precompute_elapsed = toc;
fprintf('State-transition precompute: %.2fs (%s)\n', precompute_elapsed, device_str_precompute);

%% GPU-accelerated value iteration

% L and next_state_idx are already gpuArrays (or plain arrays, if
% useGPU is false) from the precompute step above — no extra transfer
% needed here. next_state_idx is used purely for indexing; it was cast
% via precisionClass above, so round it back to a clean integer type
% for indexing correctness (matters mainly if precisionClass = 'single').
L_gpu = L;
next_state_idx_gpu = round(next_state_idx);
J = castfn(zeros(n_states, 1));

tol = 0.05; % delta tolerance threshold to trigger convergence
delta = inf;
iters = 0;
n_chars = 0;
y = 1; % discount

tic;
% iterative value loop to solve for J* — this is the part that actually
% benefits from the GPU: a (n_states x n_u) elementwise add + row-min,
% repeated every iteration.
while delta > tol
    [J_next, opt_control] = min(L_gpu + y*J(next_state_idx_gpu), [], 2);
    delta = gather(max(abs(J_next - J))); % gather() the scalar for the while-condition / fprintf
    J = J_next;
    iters = iters + 1;

    fprintf(repmat('\b', 1, n_chars));
    msg = sprintf('Iteration: %d  |  Delta: %.3f', iters, delta);
    n_chars = length(msg);
    fprintf(msg);
end
elapsed = toc;
if useGPU
    device_str = 'GPU';
else
    device_str = 'CPU';
end
fprintf('\nValue iteration converged in %d iterations (%.2fs, %s)\n', ...
    iters, elapsed, device_str);

% bring results back to host for plotting
J = gather(J);
opt_control = gather(opt_control);

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
% u_opt_matrix = interp2(u_opt_matrix, 2);

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
    ang_err  = wrapToPi_manual(X(:,1) - goal(1));
    vel_err  = X(:,2) - goal(2);
    is_goal  = sqrt(ang_err.^2 + vel_err.^2) < tol;
end

function C = quadratic_cost_vec(X, u, goal, tol_goal, Q, R)
    % shift theta so goal is at origin
    diff       = [wrapToPi_manual(X(:,1) - goal(1)), X(:,2) - goal(2)];
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

% GPU-safe replacements for Mapping-Toolbox wrap functions & interp1
% (elementwise arithmetic — supported for both plain doubles and gpuArray,
% no toolbox dependency, and faster than interp1/wrapToPi even on CPU)

function wrapped = wrapToPi_manual(theta)
    wrapped = mod(theta + pi, 2*pi) - pi;
end

function wrapped = wrapTo2Pi_manual(theta)
    wrapped = mod(theta, 2*pi);
end

% clamp to grid bounds then find nearest bin, assuming UNIFORM spacing
% (true for linspace-generated bins) — replaces interp1-based val_to_idx
function idx = val_to_idx_uniform(val, bins)
    val_clamped = min(max(val, bins(1)), bins(end));
    step = bins(2) - bins(1);
    idx = round((val_clamped - bins(1)) / step) + 1;
end
