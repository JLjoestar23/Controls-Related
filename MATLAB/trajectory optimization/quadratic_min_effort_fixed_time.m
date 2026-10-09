function quadratic_min_effort_fixed_time()
    dt = 0.01;
    Ad = eye(2) + dt*[0, 1; 0, 0];
    Bd = dt*[0; 1];
    x0 = [-2; 0];
    xf = [0; 0];
    nx = 2; nu = 1;

    % --- Step 1: find minimum-time N via feasibility bisection (unchanged) ---
    N_lo = 2; N_hi = 1000;
    while N_hi - N_lo > 1
        N_mid = floor((N_lo + N_hi)/2);
        if is_feasible(N_mid, Ad, Bd, x0, xf)
            N_hi = N_mid;
        else
            N_lo = N_mid;
        end
    end
    N_min = N_hi;
    fprintf('N_min = %d, t_f = %.4f s\n', N_min, (N_min-1)*dt);

    % --- Step 2: re-solve at (N_min) as a QP, minimizing effort ---
    N = N_min;
    [x_sol, u_sol] = solve_quadprog(N, Ad, Bd, x0, xf);

    t_x = (0:N-1)*dt;
    t_u = (0:N-2)*dt;

    figure;
    subplot(3,1,1);
    plot(x_sol(1,:), x_sol(2,:), 'LineWidth', 2);
    xlabel('q'); ylabel('qdot'); title('Phase Portrait'); grid on;

    subplot(3,1,2);
    plot(t_x, x_sol(1,:), t_x, x_sol(2,:), 'LineWidth', 2);
    xlabel('t (s)'); ylabel('state'); legend('q','qdot'); title('States vs Time'); grid on;

    subplot(3,1,3);
    plot(t_u, u_sol, 'LineWidth', 2);
    ylim([-1.1, 1.1]);
    xlabel('t (s)'); ylabel('u'); title('Control Input (min effort)'); grid on;
end

function feasible = is_feasible(N, A, B, x0, xf)
    [Aeq, beq, lb, ub, nVars] = build_constraints(N, A, B, x0, xf);
    opts = optimoptions('linprog', 'Display', 'none');
    [~, ~, exitflag] = linprog(zeros(nVars,1), [], [], Aeq, beq, lb, ub, opts);
    feasible = (exitflag == 1);
end

function [x_sol, u_sol] = solve_quadprog(N, A, B, x0, xf)
    [Aeq, beq, lb, ub, nVars, nX, nu] = build_constraints(N, A, B, x0, xf);

    H = zeros(nVars);
    H(nX+1:end, nX+1:end) = 2*eye(nVars - nX);   % quadprog uses (1/2)x'Hx, so 2*I gives sum(u.^2)
    fcost = zeros(nVars,1);

    opts = optimoptions('quadprog', 'Display', 'none');
    z = quadprog(H, fcost, [], [], Aeq, beq, lb, ub, [], opts);

    nx = 2;
    x_sol = reshape(z(1:nX), nx, N);
    u_sol = reshape(z(nX+1:end), nu, N-1);
end

function [Aeq, beq, lb, ub, nVars, nX, nu] = build_constraints(N, A, B, x0, xf)
    nx = 2; nu = 1;
    nX = nx*N; nU = nu*(N-1);
    nVars = nX + nU;

    Aeq = zeros(2*nx + nx*(N-1), nVars);
    beq = zeros(size(Aeq,1), 1);
    Aeq(1:nx, 1:nx) = eye(nx);
    beq(1:nx) = x0;
    idxN = (N-1)*nx + (1:nx);
    Aeq(nx+1:2*nx, idxN) = eye(nx);
    beq(nx+1:2*nx) = xf;

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
    lb(nX+1:end) = -1; ub(nX+1:end) = 1;
end