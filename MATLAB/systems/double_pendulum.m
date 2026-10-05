% double pendulum rate function
function dXdt = double_pendulum(t, X, system_params)

    m1 = system_params.m1;
    m2 = system_params.m2;
    l1 = system_params.l1;
    l2 = system_params.l2;
    g = system_params.g;

    phi1  = X(1);
    phi2  = X(2);
    dphi1 = X(3);
    dphi2 = X(4);
    q = [phi1; phi2];
    dq = [dphi1; dphi2];

    M = [(m1+m2)*l1^2,              m2*l1*l2*cos(phi2-phi1);
          m2*l1*l2*cos(phi2-phi1),  m2*l2^2];

    C = [0,                               -m2*l1*l2*dphi2*sin(phi2-phi1);
          m2*l1*l2*dphi1*sin(phi2-phi1),   0];

    gvec = -g * [(m1+m2)*l1*sin(phi1);
                m2*l2*sin(phi2)];

    B = [1 0 ; 0 1];
    
    % command acceleration to mimic single pendulum
    % ddq_des = [g*sin(phi1)/(l1+l2); g*sin(phi2)/(l1+l2)];
    % ddq_des = [pi; 0];
    
    % command acceleration based on position feedback
    % q_des = [sin(t); 0.5*sin(t)];
    q_des = [pi; sin(t)];
    Kp = 50;
    Kd = 10;
    ddq_des = -Kp*(q-q_des) - Kd*dq;
    
    % choose whether there is input or not
    u = C*dq + gvec + M * ddq_des;
    % u = 0;
    
    accels = M \ (-C*[dphi1; dphi2] + gvec + B*u);
    dXdt   = [dphi1; dphi2; accels(1); accels(2)];
end