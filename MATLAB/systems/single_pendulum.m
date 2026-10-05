% single pendulum rate function
function dXdt = single_pendulum(t, X, system_params)

    m = system_params.m;
    l = system_params.l;
    g = 9.81;

    theta  = X(1);
    dtheta  = X(2);

    M = m*l^2;
    
    F_g = m*g*l*sin(theta);
    
    b = 0;

    % use standard PD after swingup
    if abs(pi - theta) < 0.2 && abs(dtheta) < 2.0
        Kp = 50;
        Kd = 10;
        u  = Kp*(pi - theta) - Kd*dtheta;
        u  = max(min(u, 6), -6);
    else
    %     energy shaping swingup
        E = 0.5*m*l^2*dtheta^2 - m*g*l*cos(theta);
        E_des = m*g*l;
        E_err = E - E_des;
        k = 10;
        u = -k*dtheta*E_err;
        u = max(min(u, 6), -6);
    end

    ddtheta = (-(F_g + b*dtheta) + u) / M;

    dXdt   = [dtheta; ddtheta];
end