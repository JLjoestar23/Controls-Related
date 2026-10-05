function dXdt = cartpole(t, X, system_params)
    
    mc = system_params.mc;
    mp = system_params.mp;
    l = system_params.l;
    g = 9.81;
    K = system_params.K;
    
    x = X(1);
    theta  = X(2);
    dx = X(3);
    dtheta  = X(4);

    q = [x; theta];
    dq = [dx; dtheta];

    M = [mc + mp,           mp*l*cos(theta); 
         mp*l*cos(theta),   mp*l^2];
    
    C = [0,     -mp*l*dtheta*sin(theta);
         0,     0];

    tau = [0; -mp*g*l*sin(theta)];

    B = [1; 0];
    
    if abs(wrapToPi(pi-theta)) < pi/8
        if theta <= 0
            theta_des = -pi;
        elseif theta > 0
            theta_des = pi;
        end
        u = K*([0; theta_des; 0; 0] - X); % LQR feedback
        u = max(min(u, 8), -8);
        % disp('LQR')
        % disp(t)
    else
        E = 0.5*mp*l^2*dtheta^2 - mp*g*l*cos(theta);
        E_des = 1.25*mp*g*l;
        E_err = E - E_des;
        k = 8;
        kp = 6;
        kd = 8;
        % swingup controller + PD controller for position
        u = k*dtheta*cos(theta)*E_err - kp*x - kd*dx;
        u = max(min(u, 8), -8);
    end
    
    ddq = M \ (tau + B*u - C*dq);

    dXdt = [dq; ddq];
end