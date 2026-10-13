function [X_sol, U_sol] = solve_quad_prog(N, A, B, x0, xf)

    nX = size(A, 1);
    nU = size(B, 2);
    nD = nX*N+nU*(N-1);
    nG = nX+nU;
    
    [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf);
    
    [H, f] = build_cost(N, nX, nU, 'quadratic');

    opts = optimoptions('quadprog', 'Display', 'none');
    z = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);
    
    initial_reshape = reshape([z; zeros(nU, 1)], nG, N);
    X_sol = initial_reshape(1:nX, :);
    U_sol = initial_reshape(nX+1:end, 1:N-1);
end