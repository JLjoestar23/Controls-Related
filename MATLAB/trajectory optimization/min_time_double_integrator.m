function N_min = min_time_double_integrator()
    dt = 0.005;
    A = eye(2) + dt*[0 1; 0 0];
    B = dt*[0; 1];
    x0 = [-2; 0];  xf = [0; 0];

    N_lo = 2; N_hi = 1000;   % N_hi must be a known-feasible upper bound
    while N_hi - N_lo > 1
        N_mid = floor((N_lo + N_hi)/2);
        if is_feasible(N_mid, A, B, x0, xf)
            N_hi = N_mid;
        else
            N_lo = N_mid;
        end
    end
    N_min = N_hi;
    fprintf('N_min = %d, t_f = %.4f s\n', N_min, (N_min-1)*dt);
end

function feasible = is_feasible(N, A, B, x0, xf)
    nx = 2; nu = 1;
    nX = nx*N; nU = nu*(N-1);
    nVars = nX + nU;

    Aeq = zeros(2*nx + nx*(N-1), nVars);
    beq = zeros(size(Aeq,1), 1);

    Aeq(1:nx, 1:nx) = eye(nx);              beq(1:nx) = x0;
    idxN = (N-1)*nx + (1:nx);
    Aeq(nx+1:2*nx, idxN) = eye(nx);         beq(nx+1:2*nx) = xf;

    r = 2*nx;
    for k = 1:N-1
        rows = r + (k-1)*nx + (1:nx);
        colk  = (k-1)*nx + (1:nx);
        colk1 = k*nx + (1:nx);
        colu  = nX + (k-1)*nu + (1:nu);
        Aeq(rows, colk1) =  eye(nx);
        Aeq(rows, colk)  = -A;
        Aeq(rows, colu)  = -B;
    end

    lb = -inf(nVars,1); ub = inf(nVars,1);
    lb(nX+1:end) = -1;  ub(nX+1:end) = 1;    % |u| <= 1

    f = zeros(nVars,1);                      % feasibility only — no objective
    opts = optimoptions('linprog','Display','none');
    [~,~,exitflag] = linprog(f, [], [], Aeq, beq, lb, ub, opts);
    feasible = (exitflag == 1);
end