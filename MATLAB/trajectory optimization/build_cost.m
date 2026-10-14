function [H, f] = build_cost(N, nX, nU, cost_type)
    nD = nX*N+nU*(N-1);
    if strcmpi(cost_type, 'quadratic')
        H = 2*eye(nD);
        f = zeros(nD, 1);
    elseif strcmpi(cost_type, 'linear')
        H = zeros(nD);
        f = zeros(nD, 1);
    elseif strcmpi(cost_type, 'mineffort')
        nz = nX + nU;
        H = sparse(nD, nD);
        for k = 1:N
           x_idx = (k-1)*nz + (1:nX);
           H(x_idx, x_idx) = 1e-6*speye(nX);
        end
        for k = 1:N-1
           u_idx = (k-1)*nz + nX + (1:nU);
           H(u_idx, u_idx) = 2*speye(nU);
        end
       f = zeros(nD, 1);
    end
end