function [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf)
    % N; % number of discrete timesteps
    % nX; % state dimensions
    % nU; % input dimensions
    
    % A; % discrete-time state matrix
    % B; % discrete-time input matrix
    % x0; % system ICs

    nD = nX*N + nU*(N-1); % total number of decision variables
    nG = nX+nU; % number of decision variables per "group"
    
    % initialize equality constraint matrix and vector
    Aeq = zeros(nG*N, nD);
    beq = zeros(nG*N, 1);
    
    % setup for system initial conditions
    Aeq(1:nX, 1:nX) = eye(2);
    beq(1:nX) = x0;
    
    % dynamics equality constraints: x[k+1] = Ax[k] + Bu[k]
    for k=1:N-1
        row = nX + (k-1)*nX + (1:nX);
        x_col = (k-1)*nG + (1:nX);
        x_next_col = k*nG + (1:nX);
        u_col = x_col(end) + (1:nU);
        Aeq(row, x_col) = -A;
        Aeq(row, x_next_col) = eye(2);
        Aeq(row, u_col) = -B;
    end

    % setup for system terminal conditions
    Aeq(end-(nX-1):end, end-(nX-1):end) = eye(2);
    beq(end-(nX-1):end) = xf;
    
    % setup lower and upper bounds for states and systems
    lb = -inf(nD,1);
    ub = inf(nD,1);

    % define bounded input
    for k = 1:N-1
        u_idx = (k-1)*(nX+nU) + nX + (1:nU);
        lb(u_idx) = -1;  ub(u_idx) = 1;
    end
end