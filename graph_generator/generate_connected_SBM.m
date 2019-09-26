% Generate a connected graph from stochastic block model (SBM)
% Input:  #nodes = N_v; #modules = N_m; edge prob. across/within modules p1/p2
% Output: adjacency matrix A, combinatorial Laplacian L
% Yu Zhu, Rice ECE, 12/04/2018
function [A,L] = generate_connected_SBM(N_v,N_m,p1,p2)
    n = N_v/N_m; % Number of nodes in each module
    flag_connected = 0;
    while flag_connected ==0
        A = rand(N_v) < p1;     
        A = A.*(1-kron(eye(N_m),ones(n))); % Set diagonal blocks zero
        D = rand(n) < p2;     
        for i = 1:(N_m-1)
            temp = rand(n) < p2; 
            D = blkdiag(D,temp);
        end
        A = A + D;
        A = triu(A,1);
        A = A + A';
        L = diag(sum(A)) - A;
        if sum(abs(eig(L)) <= 1.00e-06) == 1
            flag_connected = 1; % Check whether the graph is connected
        end
    end
end


