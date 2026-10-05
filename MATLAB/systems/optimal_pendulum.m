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
X0    = [0.01; 0];
opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',0.01);
[t, X] = ode45(@(t,X) single_pendulum(t, X, system_params, opt_params), tspan, X0, opts);

%%

figure;
plot(t, X);
grid on;
legend('theta', 'dtheta');

single_pendulum_anim(X, t, tspan, system_params)

%%
function dXdt = single_pendulum(t, X, system_params, opt_params)
    m  = system_params.m;
    l  = system_params.l;
    g  = system_params.g;

    u_opt_mat   = opt_params.u_opt_mat;
    theta_bins  = opt_params.theta_bins;
    dtheta_bins = opt_params.dtheta_bins;
    d_theta     = opt_params.d_theta;
    d_dtheta    = opt_params.d_dtheta;

    theta         = X(1);
    dtheta        = X(2);
    theta_wrapped = wrapTo2Pi(theta);

    if abs(pi - theta_wrapped) < 0.05
        Kp = 50; Kd = 10;
        u  = Kp*(pi - theta_wrapped) - Kd*dtheta;
        u  = max(min(u, 6), -6);
    else
        th_c  = min(max(theta_wrapped, theta_bins(1)),  theta_bins(end));
        dth_c = min(max(dtheta,        dtheta_bins(1)), dtheta_bins(end));

        ti  = min(max(round((th_c  - theta_bins(1))  / d_theta)  + 1, 1), numel(theta_bins));
        dti = min(max(round((dth_c - dtheta_bins(1)) / d_dtheta) + 1, 1), numel(dtheta_bins));

        u = u_opt_mat(ti, dti);
    end

    ddtheta = (-m*g*l*sin(theta) + u) / (m*l^2);
    dXdt    = [dtheta; ddtheta];
end