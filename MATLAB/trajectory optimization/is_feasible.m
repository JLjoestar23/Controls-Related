function feasible = is_feasible(N, A, B, x0, xf)
    nX = size(A, 1);
    nU = size(B, 2);
    nD = nX*N+nU*(N-1);
    % nG = nX+nU;
    
    [Aeq, beq, lb, ub] = build_constraints(N, nX, nU, A, B, x0, xf);

    opts = optimoptions('linprog', 'Display', 'none');
    [~, ~, exit_flag] = linprog(zeros(nD,1), [], [], Aeq, beq, lb, ub, opts);
    feasible = (exit_flag == 1);
end