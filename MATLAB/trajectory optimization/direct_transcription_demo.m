%% Direct Transcription Demo for a Double Integrator System
% Demonstrates the implementation and solving for the optimal trajectory
% and control policy for a arbitrary double integrator system

%% Establish System Dynamics

dt = 0.01;
Ad = eye(2) + dt * [0, 1; 0, 0];
Bd = dt * [0; 1];

function [result, prog, X, u] = solve_for_fixed_horizon(N)
    
    % decision variable structure:
    % N state variables and N - 1 input variables
    
    % define minimum time cost function
    f_cost = zeros(1, 2*N-1);
    f_cost(1:N) = 1;

    for k=2:N-1
        Aeq(:, k+1) = Ad

    linprog(f, Aeq, beq)
end

N = 284;