function [H, f] = build_cost(N, nX, nU, cost_type)
    if strcmpi(cost_type, 'quadratic')
        nD = nX*N+nU*(N-1);
        H = 2*eye(nD);
        f = zeros(nD, 1);
    end
end