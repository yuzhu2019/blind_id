% Generate a connected ER graph
% Input:  graph size N, edge prob. p
% Output: adjacency matrix A, combinatorial Laplacian L
% Yu Zhu, Rice ECE, 12/04/2018
function [A,L] = generate_connected_ER(N,p)
    flag_connected = 0;
    while flag_connected ==0
        A = rand(N) < p;
        A = triu(A,1);
        A = A + A'; 
        L = diag(sum(A)) - A; 
        if sum(abs(eig(L)) <= 1.00e-06) == 1
            flag_connected = 1; % Check whether the graph is connected
        end
    end
end
