function dXdt = single_pendulum_value_iteration(t, X, system_params, opt_params)
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
        th_c  = min(max(theta_wrapped, theta_bins(1)), theta_bins(end));
        dth_c = min(max(dtheta, dtheta_bins(1)), dtheta_bins(end));

        ti  = min(max(round((th_c - theta_bins(1)) / d_theta) + 1, 1), numel(theta_bins));
        dti = min(max(round((dth_c - dtheta_bins(1)) / d_dtheta) + 1, 1), numel(dtheta_bins));
        
        u = u_opt_mat(ti, dti);
    end

    ddtheta = (-m*g*l*sin(theta) + u) / (m*l^2);
    dXdt    = [dtheta; ddtheta];
end