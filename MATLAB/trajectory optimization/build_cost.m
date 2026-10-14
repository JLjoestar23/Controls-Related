function [H, f] = build_cost(N, nX, nU, cost_type)
    nD = nX*N+nU*(N-1);
    if strcmpi(cost_type, 'quadratic')
        H = 2*eye(nD);
        f = zeros(nD, 1);
    elseif strcmpi(cost_type, 'linear')
        H = zeros(nD);
        f = zeros(nD, 1);
    end
end