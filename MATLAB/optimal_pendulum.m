%% Single Pendulum w/ Optimal Swing-up
% Parameters
system_params.m = 2; 
system_params.l = 1;
system_params.g = 9.81;

opt_params.u_opt_mat = u_opt_matrix;
opt_params.theta_bins = theta_bins;
opt_params.dtheta_bins = dtheta_bins;
opt_params.d_theta  = theta_bins(2)  - theta_bins(1);
opt_params.d_dtheta = dtheta_bins(2) - dtheta_bins(1);

% Simulate
tspan = [0 15];
X0    = [1e-3; 0];
opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',0.01);
[t, X] = ode45(@(t,X) single_pendulum_value_iteration(t, X, system_params, opt_params), tspan, X0, opts);

%%

figure;
plot(t, X);
grid on;
legend('theta', 'dtheta');

single_pendulum_anim(X, t, tspan, system_params)